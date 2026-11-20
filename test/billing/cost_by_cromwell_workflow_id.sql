-- Cost of specific Cromwell workflows, for runs from before the CI labels were added.
-- Workflow IDs are in the CI log ("[<Workflow>/<test>]     Workflow id: <uuid>");
-- in billing they appear as the label cromwell-workflow-id=cromwell-<uuid>.
-- Replace <TABLE> with broad-gp-hydrogen.billing_export.gcp_billing_export_resource_v1_<ACCOUNT_ID>
-- (billing account ID with "-" replaced by "_"). Net cost = cost + credits (e.g. spot/sustained-use discounts).
-- The export is partitioned by _PARTITIONTIME; the extra 2 days cover late-arriving rows.
SELECT
  (SELECT value FROM UNNEST(labels) WHERE key = 'cromwell-workflow-id') AS workflow_id,
  ROUND(SUM(cost + IFNULL((SELECT SUM(c.amount) FROM UNNEST(credits) c), 0)), 2) AS net_cost
FROM `<TABLE>`
WHERE DATE(_PARTITIONTIME) >= DATE_SUB(CURRENT_DATE(), INTERVAL 14 DAY)
  AND project.id = 'broad-gp-hydrogen'
  AND EXISTS (SELECT 1 FROM UNNEST(labels) l
              WHERE l.key = 'cromwell-workflow-id'
                AND l.value IN ('cromwell-<uuid1>', 'cromwell-<uuid2>'))
GROUP BY workflow_id;
