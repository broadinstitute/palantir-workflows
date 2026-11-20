-- Cost of each CI run (label ci-run=run-<GitHub run id>), split into compute and disk.
-- Bucket storage has no ci-run label, so it is not included; see daily_ci_cost_by_service.sql.
-- GitHub run: https://github.com/broadinstitute/palantir-workflows/actions/runs/<ci_run without "run-">
-- TODO: check the sku.description patterns against real rows.
-- Replace <TABLE> with broad-gp-hydrogen.billing_export.gcp_billing_export_resource_v1_<ACCOUNT_ID>
-- (billing account ID with "-" replaced by "_"). Net cost = cost + credits (e.g. spot/sustained-use discounts).
-- The export is partitioned by _PARTITIONTIME; the extra 2 days cover late-arriving rows.
DECLARE days INT64 DEFAULT 30;
WITH rows AS (
  SELECT
    (SELECT value FROM UNNEST(labels) WHERE key = 'ci-run') AS ci_run,
    (SELECT value FROM UNNEST(labels) WHERE key = 'ci-trigger') AS ci_trigger,
    usage_start_time,
    CASE
      WHEN REGEXP_CONTAINS(sku.description, r'(?i)\bPD\b|Persistent Disk|Storage PD') THEN 'disk'
      WHEN service.description = 'Compute Engine' THEN 'compute'
      ELSE 'other'
    END AS category,
    cost + IFNULL((SELECT SUM(c.amount) FROM UNNEST(credits) c), 0) AS net_cost
  FROM `<TABLE>`
  WHERE DATE(_PARTITIONTIME) >= DATE_SUB(CURRENT_DATE(), INTERVAL days + 2 DAY)
    AND project.id = 'broad-gp-hydrogen'
    AND EXISTS (SELECT 1 FROM UNNEST(labels) l WHERE l.key = 'palantir-ci' AND l.value = 'true')
)
SELECT
  ci_run, ci_trigger,
  MIN(usage_start_time) AS started,
  ROUND(SUM(IF(category = 'compute', net_cost, 0)), 2) AS compute,
  ROUND(SUM(IF(category = 'disk', net_cost, 0)), 2) AS disk,
  ROUND(SUM(IF(category = 'other', net_cost, 0)), 2) AS other,
  ROUND(SUM(net_cost), 2) AS total
FROM rows
GROUP BY ci_run, ci_trigger
ORDER BY started DESC;
