#!/usr/bin/env bash
#
# Télécharge les jeux du manifeste et vérifie leur empreinte.
#
# Deux règles, non négociables :
#   - une ligne `unresolved` fait échouer le fetch. On ne devine jamais une accession.
#   - une empreinte inscrite qui diverge fait échouer le fetch. On ne réécrit jamais
#     silencieusement un sha256 : il faut passer par `lock.sh`, qui est explicite.
#
# Variables : MANIFEST, RAW, TIERS (liste séparée par des espaces).

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MANIFEST="${MANIFEST:-$ROOT/data/manifest.tsv}"
RAW="${RAW:-$ROOT/data/raw}"
TIERS="${TIERS:-core}"

[[ -f $MANIFEST ]] || { echo "manifeste introuvable : $MANIFEST" >&2; exit 2; }

want_tier() {
  local t=$1 w
  for w in $TIERS; do [[ $w == all || $w == "$t" ]] && return 0; done
  return 1
}

n_ok=0 n_new=0 n_skip=0 n_bad=0 n_unresolved=0

while IFS=$'\t' read -r id tier status license url sha bytes notes || [[ -n ${id:-} ]]; do
  [[ -z ${id:-} || $id == \#* || $id == id ]] && continue
  want_tier "$tier" || { n_skip=$((n_skip + 1)); continue; }

  if [[ $status == unresolved || $url == UNRESOLVED ]]; then
    printf '  UNRESOLVED  %-34s %s\n' "$id" "$notes" >&2
    n_unresolved=$((n_unresolved + 1))
    continue
  fi

  dest="$RAW/$tier/$id"
  mkdir -p "$(dirname "$dest")"

  if [[ -f $dest ]]; then
    have=$(sha256sum "$dest" | cut -d' ' -f1)
    if [[ $sha == PENDING ]]; then
      printf '  PENDING     %-34s %s\n' "$id" "$have"
      n_new=$((n_new + 1))
    elif [[ $have == "$sha" ]]; then
      printf '  ok          %-34s\n' "$id"
      n_ok=$((n_ok + 1))
    else
      printf '  MISMATCH    %-34s attendu %s, obtenu %s\n' "$id" "$sha" "$have" >&2
      n_bad=$((n_bad + 1))
    fi
    continue
  fi

  printf '  fetch       %-34s %s\n' "$id" "$url"
  tmp="$dest.part"
  if ! curl -fsSL --retry 4 --retry-delay 2 --retry-connrefused -o "$tmp" "$url"; then
    printf '  FAILED      %-34s téléchargement impossible\n' "$id" >&2
    rm -f "$tmp"
    n_bad=$((n_bad + 1))
    continue
  fi

  have=$(sha256sum "$tmp" | cut -d' ' -f1)
  if [[ $sha != PENDING && $have != "$sha" ]]; then
    printf '  MISMATCH    %-34s attendu %s, obtenu %s\n' "$id" "$sha" "$have" >&2
    rm -f "$tmp"
    n_bad=$((n_bad + 1))
    continue
  fi
  mv "$tmp" "$dest"

  if [[ $sha == PENDING ]]; then
    printf '  PENDING     %-34s %s\n' "$id" "$have"
    n_new=$((n_new + 1))
  else
    printf '  ok          %-34s\n' "$id"
    n_ok=$((n_ok + 1))
  fi
done < "$MANIFEST"

echo
printf 'vérifiés %d · à verrouiller %d · hors tier %d · non résolus %d · en échec %d\n' \
  "$n_ok" "$n_new" "$n_skip" "$n_unresolved" "$n_bad"

if ((n_unresolved > 0)); then
  echo "→ figer les accessions manquantes dans le manifeste avant d'aller plus loin." >&2
fi
if ((n_new > 0)); then
  echo "→ 'make data-lock' pour inscrire les empreintes des fichiers nouvellement récupérés." >&2
fi

((n_bad == 0 && n_unresolved == 0))
