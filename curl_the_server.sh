curl -s http://Jesse:8787/status | jq -r '
  "\(.title) — \(.subtitle)",
  "Stage: \(.stage)",
  "",
  (.sections[] |
    "== \(.title) ==",
    (.metrics[] | "  \(.label + ":") \(.value)")
  )
'
