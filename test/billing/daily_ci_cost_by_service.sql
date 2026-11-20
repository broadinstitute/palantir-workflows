-- Daily cost of the WDL CI tests (label palantir-ci=true), by service, over the last `days` days.
-- Includes the test buckets once they carry the palantir-ci label.
-- Replace <TABLE> with broad-gp-hydrogen.billing_export.gcp_billing_export_resource_v1_<ACCOUNT_ID>
-- (billing account ID with "-" replaced by "_"). Net cost = cost + credits (e.g. spot/sustained-use discounts).
-- The export is partitioned by _PARTITIONTIME; the extra 2 days cover late-arriving rows.
DECLARE days INT64 DEFAULT 30;
SELECT
  DATE(usage_start_time) AS day,
  service.description AS service,
  ROUND(SUM(cost + IFNULL((SELECT SUM(c.amount) FROM UNNEST(credits) c), 0)), 2) AS net_cost,
  ANY_VALUE(currency) AS currency
FROM `<TABLE>`
WHERE DATE(_PARTITIONTIME) >= DATE_SUB(CURRENT_DATE(), INTERVAL days + 2 DAY)
  AND usage_start_time >= TIMESTAMP(DATE_SUB(CURRENT_DATE(), INTERVAL days DAY))
  AND project.id = 'broad-gp-hydrogen'
  AND EXISTS (SELECT 1 FROM UNNEST(labels) l WHERE l.key = 'palantir-ci' AND l.value = 'true')
GROUP BY day, service
ORDER BY day DESC, net_cost DESC;
