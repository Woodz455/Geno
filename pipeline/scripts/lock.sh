#!/usr/bin/env bash
#
# Inscrit dans le manifeste l'empreinte et la taille des fichiers déjà récupérés,
# et bascule leur statut en `ok`.
#
# Étape séparée et explicite : un fetch ne réécrit jamais un sha256 de lui-même.
# Verrouiller, c'est affirmer « ce fichier est celui que je veux ». Ça se fait à la main.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MANIFEST="${MANIFEST:-$ROOT/data/manifest.tsv}"
RAW="${RAW:-$ROOT/data/raw}"

[[ -f $MANIFEST ]] || { echo "manifeste introuvable : $MANIFEST" >&2; exit 2; }

tmp=$(mktemp)
locked=0

while IFS= read -r line || [[ -n $line ]]; do
  if [[ -z $line || $line == \#* || $line == id$'\t'* ]]; then
    printf '%s\n' "$line" >> "$tmp"
    continue
  fi

  IFS=$'\t' read -r id tier status license url sha bytes notes <<< "$line"
  dest="$RAW/$tier/$id"

  if [[ $sha == PENDING && -f $dest ]]; then
    new_sha=$(sha256sum "$dest" | cut -d' ' -f1)
    new_bytes=$(wc -c < "$dest" | tr -d ' ')
    printf '%s\t%s\tok\t%s\t%s\t%s\t%s\t%s\n' \
      "$id" "$tier" "$license" "$url" "$new_sha" "$new_bytes" "$notes" >> "$tmp"
    printf '  verrouillé  %-34s %s  (%s octets)\n' "$id" "$new_sha" "$new_bytes"
    locked=$((locked + 1))
  else
    printf '%s\n' "$line" >> "$tmp"
  fi
done < "$MANIFEST"

mv "$tmp" "$MANIFEST"
echo
printf '%d entrée(s) verrouillée(s).\n' "$locked"
