# Delivering data from ICA to S3

`deliver_to_s3.py` is a small helper script that copies a folder's worth of
files from Illumina Connected Analytics (ICA) into
`s3://cell-programs-data/screens/...`, where the data is kept long-term. You
run it from a terminal, answer a few questions, and it does the rest. You do
not need to know anything about AWS or S3 to use it — just ICA.

(If you want to understand *how* the delivery pipeline works behind the
scenes — the S3 buckets, replication, permissions, etc. — see
[`ica-s3-delivery-setup.md`](ica-s3-delivery-setup.md) in this same folder.
You don't need to read that to just deliver some data.)

## 1. One-time setup

### 1.1 Get the script onto your computer

If you already have this repository checked out, you already have the
script — it's the `deliver_to_s3.py` file in this same folder. Otherwise, ask
whoever shared this with you for the `ica_tools/delivery/` folder (both
`deliver_to_s3.py` and this `README.md`).

### 1.2 Make sure you have Python

The script needs Python 3.8 or newer, and nothing else — no extra packages
to install. Most Macs and Linux machines already have it. To check, open a
terminal and run:

```
python3 --version
```

If that prints a version number (3.8 or higher), you're set. If it says
"command not found" or shows an older version, ask an administrator for
help installing Python 3.

### 1.3 Get an ICA API key

The script needs to act on your behalf in ICA, which requires an API key —
a long password-like string tied to your ICA account.

1. Log in to [ICA](https://ica.illumina.com).
2. Click your user icon (top right) → **Account Settings** (or similar,
   depending on your ICA version) → look for **API Keys**.
3. Create a new API key and copy it somewhere safe for the next step.

**Treat this key like a password.** Anyone with it can act as you in ICA.
Don't share it, paste it into chat messages, or commit it to a code
repository.

### 1.4 Save the API key where the script can find it

Save the key as plain text, and nothing else (no quotes, no extra text), in
**one** of these two locations:

- A file named `illumina-api-key.txt` in this same folder
  (`ica_tools/delivery/illumina-api-key.txt`), **or**
- A file at `~/.icav2/api_key.txt` (in your home folder, inside a folder
  called `.icav2`) — this location is shared with other ICA tools in this
  repository, so if you've already set one up for those, you don't need to
  do anything else.

To create the first option from a terminal, `cd` into this folder and run:

```
echo "YOUR_API_KEY_HERE" > illumina-api-key.txt
```

replacing `YOUR_API_KEY_HERE` with the key you copied. If this file lives
inside a checked-out copy of this repository, it will **not** accidentally
get committed to git — it's already excluded.

## 2. Running the script

Open a terminal, navigate to this folder, and run:

```
python3 deliver_to_s3.py
```

(On some systems the command is `python` instead of `python3`.)

The script will ask you a few questions:

### 2.1 Which ICA folder to deliver

It asks for the **ICA folder ID** of the folder you want to deliver — not
its name, but an ID that looks like `fol.5e5f...` (starts with `fol.`
followed by a string of letters and numbers).

To find it in ICA: navigate to the folder, then look at its **Data
details** panel — the ID is shown there. Copy it exactly.

The script delivers the **contents** of that folder (all the files and
subfolders inside it), not the folder itself.

### 2.2 Where it should go

Next, it asks for a destination path under `screens/`, for example:

```
PotC_Screen_1/sample_1/processing_output
```

This becomes the folder structure your files will have in S3, i.e. they'll
end up at:

```
s3://cell-programs-data/screens/PotC_Screen_1/sample_1/processing_output/
```

Pick a destination that doesn't already exist, so your delivery doesn't mix
with someone else's — using a unique combination such as
screen name / sample name / a label for what the data is (e.g.
`processing_output`) usually works well. You don't need to create these
folders beforehand; the script creates them for you if they don't exist.

Only letters, digits, dots (`.`), dashes (`-`), and underscores (`_`) are
allowed in each part of the path (no spaces or other punctuation).

### 2.3 Confirm before anything happens

Before copying anything, the script shows you a summary: where the files
are coming from, where they're going, how many files/folders are involved,
and any warnings (for example, if the destination already has files in it,
or if some source files aren't available to copy). **Nothing is copied
until you type `yes`** at this point — if anything looks wrong, just answer
something else (or press Ctrl+C) to cancel, and nothing will happen.

### 2.4 What happens after you confirm

Once you type `yes`, the script tells ICA to start copying your files. This
copy happens inside ICA/AWS, in the background — it does **not** depend on
your computer staying on or connected. The script prints:

- a batch job ID (or IDs, if there were a lot of files),
- a link you can open in your browser to watch the copy's progress in ICA
  (under that project's **Activity → Batch Jobs**), and
- the final S3 location where your files will show up once the copy
  finishes (usually within a few minutes, longer for very large deliveries).

You can close the terminal immediately after this — the copy keeps running
on ICA's side regardless.

## 3. Troubleshooting

- **"API key file not found"** — you haven't completed step 1.4 above, or
  saved the file in the wrong place/name. Double check the exact file name
  and location.
- **"ICA rejected the API key (401)"** — the key is wrong, expired, or was
  revoked. Create a new one (step 1.3) and save it again (step 1.4).
- **"Folder ... was not found in any ICA project you have access to"** —
  double-check the folder ID, and make sure your ICA account actually has
  access to the project that folder is in.
- **It says some files were skipped** — files that aren't in an "available"
  state in ICA (for example, archived files) are skipped automatically; the
  summary screen tells you which ones. Files that already exist at the
  destination are also skipped, rather than overwritten — this makes it
  safe to run the script again with the same destination if a previous
  delivery didn't fully finish.
- **Running it again** — it's always safe to run the script again with the
  same source and destination. It will only copy whatever hasn't already
  been delivered there.
- **Anything else** — note down the exact error message and ask an
  administrator for help.
