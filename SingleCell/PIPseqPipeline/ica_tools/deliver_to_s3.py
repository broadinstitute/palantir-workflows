#!/usr/bin/env python3
"""
Deliver an ICA folder to s3://cell-programs-data/screens/...

Copies the CONTENTS of one ICA folder (from any project you can access) into a
folder of the ICA delivery project. Everything in the delivery project is
replicated automatically to s3://cell-programs-data/screens/.

Example: the contents of samples/sample_id_123/ copied to destination
"PotC_Screen_1/sample_1/processing_output" end up in
s3://cell-programs-data/screens/PotC_Screen_1/sample_1/processing_output/

Requirements
  * Python 3.8 or newer (no extra packages needed).
  * Your ICA API key, saved on its own in either "illumina-api-key.txt" in the
    same folder as this script, or "~/.icav2/api_key.txt". Keep that file private.

Run it with:  python deliver_to_s3.py
"""

import json
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

try:
    # Importing this hooks into input() to support arrow keys/history editing at the
    # prompts below. Not available on Windows, where input() just works without it.
    import readline  # noqa: F401
except ImportError:
    pass

# --------------------------------------------------------------------------
# Configuration (set once by an administrator)
# --------------------------------------------------------------------------
ICA_BASE_URL = "https://ica.illumina.com/ica/rest"
# ID of the ICA delivery project (Project > Details > URN, the part after "project:").
DELIVERY_PROJECT_ID = "159ebcf0-88b6-4dc3-ab99-e2378001ea32"
# Where the delivery project's root ends up in S3 (shown to the user only).
FINAL_S3_PREFIX = "s3://cell-programs-data/screens/"
# Checked in order; the first one that exists is used.
API_KEY_FILES = [
    Path(__file__).resolve().parent / "illumina-api-key.txt",
    Path("~/.icav2/api_key.txt").expanduser(),
]

COPY_BATCH_SIZE = 100  # number of items per copy request
MEDIA_TYPE = "application/vnd.illumina.v3+json"

FOLDER_ID_RE = re.compile(r"^fol\.[0-9a-f]{32}$")
# Allowed folder names: letters, digits, dot, dash, underscore (safe in S3 keys).
SEGMENT_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


class DeliveryError(Exception):
    """An error with a message meant for the user."""


# --------------------------------------------------------------------------
# ICA API access
# --------------------------------------------------------------------------
class IcaApi:
    def __init__(self, api_key):
        self.api_key = api_key

    def request(self, method, path, params=None, body=None, none_on=()):
        """Call the ICA API. Returns parsed JSON, or None if the HTTP status is in none_on."""
        url = ICA_BASE_URL + path
        if params:
            url += "?" + urllib.parse.urlencode(
                {k: v for k, v in params.items() if v is not None}
            )
        data = json.dumps(body).encode("utf-8") if body is not None else None
        req = urllib.request.Request(url, data=data, method=method)
        req.add_header("X-API-Key", self.api_key)
        req.add_header("Accept", f"{MEDIA_TYPE}, application/json;q=0.9, */*;q=0.5")
        if data is not None:
            req.add_header("Content-Type", MEDIA_TYPE)
        try:
            with urllib.request.urlopen(req, timeout=120) as resp:
                raw = resp.read()
                return json.loads(raw) if raw else {}
        except urllib.error.HTTPError as e:
            if e.code in none_on:
                return None
            raise DeliveryError(_describe_http_error(e)) from None
        except urllib.error.URLError as e:
            raise DeliveryError(
                f"Could not connect to ICA ({e.reason}). Check your internet connection."
            ) from None

    def paged(self, path, params):
        """Yield all items of a paged list endpoint (cursor-based pagination)."""
        params = dict(params, pageSize="1000")
        seen_tokens = set()
        while True:
            resp = self.request("GET", path, params=params)
            yield from resp.get("items", [])
            token = resp.get("nextPageToken")
            if not token or token in seen_tokens:
                return
            seen_tokens.add(token)
            params["pageToken"] = token


def _describe_http_error(e):
    detail = ""
    try:
        problem = json.loads(e.read() or b"{}")
        detail = problem.get("detail") or problem.get("title") or problem.get("message") or ""
    except (ValueError, AttributeError):
        pass
    if e.code == 401:
        return ("ICA rejected the API key (401). Check that your API key file contains "
                "a valid, non-expired key.")
    if e.code == 403:
        return f"ICA denied access (403). {detail}".strip()
    return f"ICA request failed ({e.code}). {detail}".strip()


def data_of(item):
    """Project data responses wrap the data object in 'data'; plain data responses don't."""
    return item.get("data", item) if isinstance(item, dict) else {}


def details_of(item):
    return data_of(item).get("details", {})


def norm_folder_path(path):
    return "/" + path.strip("/") + "/" if path.strip("/") else "/"


# --------------------------------------------------------------------------
# Steps
# --------------------------------------------------------------------------
def load_api_key():
    for path in API_KEY_FILES:
        if path.exists():
            key = path.read_text(encoding="utf-8").strip()
            if not key:
                raise DeliveryError(f"The API key file is empty: {path}")
            return key
    locations = "\n".join(f"  - {path}" for path in API_KEY_FILES)
    raise DeliveryError(
        "API key file not found. Save your ICA API key (and nothing else) in one "
        f"of these files:\n{locations}"
    )
    return key


def get_project_name(api, project_id):
    project = api.request("GET", f"/api/projects/{project_id}")
    return project.get("name", project_id)


def find_source_folder(api, folder_id):
    """Return (project_id, folder) for a folder ID, searching projects if needed."""
    # Fast path: look the data up directly, then confirm we can access its project.
    found = api.request("GET", f"/api/data/{folder_id}", none_on=(400, 403, 404))
    if found:
        project_id = details_of(found).get("owningProjectId")
        if project_id:
            in_project = api.request(
                "GET", f"/api/projects/{project_id}/data/{folder_id}", none_on=(400, 403, 404)
            )
            if in_project:
                return project_id, in_project

    # Fallback: check every project the user can access.
    for project in api.paged("/api/projects", {}):
        project_id = project.get("id")
        if not project_id:
            continue
        in_project = api.request(
            "GET", f"/api/projects/{project_id}/data/{folder_id}", none_on=(400, 403, 404)
        )
        if in_project:
            return project_id, in_project

    raise DeliveryError(
        f"Folder {folder_id} was not found in any ICA project you have access to."
    )


def list_children(api, project_id, folder_id):
    return list(api.paged(f"/api/projects/{project_id}/data", {"parentFolderId": folder_id}))


def find_folder_by_path(api, project_id, path):
    """Return the folder at exactly this path in the project, or None."""
    resp = api.request(
        "GET",
        f"/api/projects/{project_id}/data",
        params={
            "filePath": path,
            "filePathMatchMode": "FULL_CASE_INSENSITIVE",
            "type": "FOLDER",
            "pageSize": "100",
        },
    )
    for item in resp.get("items", []):
        found_path = norm_folder_path(details_of(item).get("path", ""))
        if found_path == path:
            return item
        if found_path.lower() == path.lower():
            raise DeliveryError(
                f"The delivery project already has a folder '{found_path}', which differs "
                f"from '{path}' only in upper/lower case. Please use the exact existing "
                "name or choose a different destination."
            )
    return None


def parse_destination(text):
    """Turn user input into a list of folder names, e.g. ['PotC_Screen_1', 'sample_1', 'processing_output']."""
    t = text.strip().strip("\"'").replace("\\", "/")
    for prefix in (FINAL_S3_PREFIX, FINAL_S3_PREFIX[len("s3://"):]):
        if t.lower().startswith(prefix.lower()):
            t = t[len(prefix):]
            break
    segments = [s for s in t.split("/") if s]
    if not segments:
        raise DeliveryError("Please enter a destination folder, e.g. PotC_Screen_1/sample_1/processing_output")
    if segments[0].lower() == "screens":
        raise DeliveryError(
            "Leave out 'screens/' at the start; it is added automatically. "
            "Example: PotC_Screen_1/sample_1/processing_output"
        )
    for seg in segments:
        if seg in (".", "..") or not SEGMENT_RE.match(seg):
            raise DeliveryError(
                f"'{seg}' is not an allowed folder name. Use only letters, digits, "
                "'.', '-' and '_' (and don't start with '.', '-' or '_')."
            )
    return segments


def ensure_destination_folder(api, project_id, segments):
    """Create the destination folders (if needed) and return the last folder's ID."""
    parent_path = "/"
    folder = None
    for seg in segments:
        path = parent_path + seg + "/"
        folder = find_folder_by_path(api, project_id, path)
        if folder is None:
            folder = api.request(
                "POST",
                f"/api/projects/{project_id}/data",
                body={"name": seg, "folderPath": parent_path, "dataType": "FOLDER"},
            )
            print(f"  Created folder {path}")
        parent_path = path
    folder_id = data_of(folder).get("id")
    if not folder_id:
        raise DeliveryError("Could not determine the ID of the destination folder.")
    return folder_id


def start_copy(api, project_id, item_ids, destination_folder_id):
    batch_ids = []
    for start in range(0, len(item_ids), COPY_BATCH_SIZE):
        chunk = item_ids[start:start + COPY_BATCH_SIZE]
        resp = api.request(
            "POST",
            f"/api/projects/{project_id}/dataCopyBatch",
            body={
                "items": [{"dataId": i} for i in chunk],
                "destinationFolderId": destination_folder_id,
                "copyUserTags": True,
                "copyTechnicalTags": True,
                "copyInstrumentInfo": True,
                "actionOnExist": "SKIP",  # never overwrite anything already delivered
            },
        )
        batch_ids.append(resp.get("id"))
    return batch_ids


def activity_url(project_id):
    return (
        f"https://ica.illumina.com/ica/projects/{project_id}/activity"
        "?tabsheet-projectactivityview=tab-projectactivityview-batchjobs"
    )


def ask(prompt):
    try:
        return input(prompt)
    except EOFError:
        raise DeliveryError("No input received.") from None


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    print("=" * 70)
    print("Deliver an ICA folder to", FINAL_S3_PREFIX)
    print("=" * 70)

    if DELIVERY_PROJECT_ID.startswith("<"):
        raise DeliveryError(
            "This script is not configured yet: set DELIVERY_PROJECT_ID at the top of the script."
        )

    api = IcaApi(load_api_key())
    delivery_name = get_project_name(api, DELIVERY_PROJECT_ID)  # also tests the API key

    # 1. Source folder
    while True:
        folder_id = ask("\nICA folder ID to deliver (looks like fol.5e5f...): ").strip().lower()
        if FOLDER_ID_RE.match(folder_id):
            break
        print("  That doesn't look like a folder ID. It starts with 'fol.' followed by "
              "32 letters/digits. You can find it in ICA under Data > folder > Data details.")
    print("  Looking up folder...")
    source_project_id, source_folder = find_source_folder(api, folder_id)
    if source_project_id == DELIVERY_PROJECT_ID:
        raise DeliveryError(
            "That folder is inside the delivery project itself. Choose a folder from the "
            "project that contains the original data."
        )
    source_details = details_of(source_folder)
    if source_details.get("dataType", "FOLDER").upper() != "FOLDER":
        raise DeliveryError(f"{folder_id} is not a folder.")
    source_project_name = get_project_name(api, source_project_id)
    source_path = norm_folder_path(source_details.get("path", ""))

    children = list_children(api, source_project_id, folder_id)
    if not children:
        raise DeliveryError(f"The folder {source_path} is empty; there is nothing to deliver.")
    n_folders = sum(1 for c in children if details_of(c).get("dataType") == "FOLDER")
    n_files = len(children) - n_folders
    not_available = [
        details_of(c).get("name", "?") for c in children
        if details_of(c).get("dataType") == "FILE"
        and details_of(c).get("status", "AVAILABLE") != "AVAILABLE"
    ]

    # 2. Destination
    while True:
        dest_input = ask("\nDestination folder under screens/ (e.g. PotC_Screen_1/sample_1/processing_output): ")
        try:
            segments = parse_destination(dest_input)
            break
        except DeliveryError as e:
            print(f"  {e}")
    dest_path = "/" + "/".join(segments) + "/"
    s3_path = FINAL_S3_PREFIX + "/".join(segments) + "/"

    existing = find_folder_by_path(api, DELIVERY_PROJECT_ID, dest_path)
    existing_count = (
        len(list_children(api, DELIVERY_PROJECT_ID, data_of(existing)["id"])) if existing else 0
    )

    # 3. Confirmation
    print("\n" + "-" * 70)
    print("Please check:")
    print(f"\n  FROM  project  {source_project_name}")
    print(f"        folder   {source_path}")
    print(f"        contents {n_files} file(s) and {n_folders} folder(s) "
          "(the folder's contents are copied, not the folder itself)")
    print(f"\n  TO    project  {delivery_name}")
    print(f"        folder   {dest_path}")
    print(f"        S3       {s3_path}")
    if not_available:
        print(f"\n  NOTE: {len(not_available)} file(s) are not in status AVAILABLE "
              f"(e.g. archived) and will be skipped: {', '.join(not_available[:5])}"
              f"{' ...' if len(not_available) > 5 else ''}")
    if existing_count:
        print(f"\n  WARNING: the destination already contains {existing_count} item(s). "
              "Files that already exist there are skipped, not overwritten.")
    print("-" * 70)

    if ask("\nType 'yes' to start the copy: ").strip().lower() not in ("yes", "y"):
        print("Cancelled. Nothing was copied.")
        return

    # 4. Copy
    print("\nPreparing destination folder...")
    dest_folder_id = ensure_destination_folder(api, DELIVERY_PROJECT_ID, segments)
    print("Starting copy...")
    batch_ids = start_copy(api, DELIVERY_PROJECT_ID, [data_of(c)["id"] for c in children],
                           dest_folder_id)

    # 5. Result
    job_ids = [b for b in batch_ids if b]
    print()
    print("Copy initiated.")
    print(f"  Batch job ID(s): {', '.join(job_ids) if job_ids else 'unknown'}")
    print(f"  Check progress at: {activity_url(DELIVERY_PROJECT_ID)}")
    print(f"The files will appear in {s3_path} once the copy finishes.")


if __name__ == "__main__":
    exit_code = 0
    try:
        main()
    except DeliveryError as err:
        print(f"\nERROR: {err}")
        exit_code = 1
    except KeyboardInterrupt:
        print("\nCancelled.")
        exit_code = 1
    if sys.stdin.isatty():
        try:
            input("\nPress Enter to close...")
        except (EOFError, KeyboardInterrupt):
            pass
    sys.exit(exit_code)
