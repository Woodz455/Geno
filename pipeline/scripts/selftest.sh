#!/usr/bin/env bash
#
# Prouve que la machinerie de vérification fait ce qu'elle prétend, sans réseau :
# le manifeste de test pointe sur des URL file://.
#
# Ce qui est vérifié :
#   1. un fetch d'une entrée neuve récupère le fichier et signale son empreinte
#   2. le verrouillage inscrit la bonne empreinte et bascule le statut en ok
#   3. un second fetch est idempotent et valide
#   4. une source altérée est REJETÉE — c'est le test qui compte vraiment
#   5. une accession non résolue fait échouer le fetch au lieu d'être devinée

set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

pass=0 fail=0
check() { # check <description> <attendu:ok|ko> <code de retour>
  if { [[ $2 == ok && $3 -eq 0 ]] || [[ $2 == ko && $3 -ne 0 ]]; }; then
    printf '  \033[32mPASS\033[0m  %s\n' "$1"; pass=$((pass + 1))
  else
    printf '  \033[31mFAIL\033[0m  %s (attendu %s, code %d)\n' "$1" "$2" "$3"; fail=$((fail + 1))
  fi
}
checkeq() { # checkeq <description> <attendu> <obtenu>
  if [[ $2 == "$3" ]]; then
    printf '  \033[32mPASS\033[0m  %s\n' "$1"; pass=$((pass + 1))
  else
    printf '  \033[31mFAIL\033[0m  %s\n        attendu %s\n        obtenu  %s\n' "$1" "$2" "$3"; fail=$((fail + 1))
  fi
}

SRC="$WORK/src"; mkdir -p "$SRC"
printf 'chr1\t248956422\nchr2\t242193529\nchr3\t198295559\n' > "$SRC/sizes.txt"
TRUE_SHA=$(sha256sum "$SRC/sizes.txt" | cut -d' ' -f1)

MAN="$WORK/manifest.tsv"
RAWD="$WORK/raw"
export MANIFEST="$MAN" RAW="$RAWD" TIERS=core

header() { printf 'id\ttier\tstatus\tlicense\turl\tsha256\tbytes\tnotes\n' > "$MAN"; }
row() { printf '%s\t%s\t%s\ttest\t%s\t%s\t0\tfixture\n' "$1" core "$2" "$3" "$4" >> "$MAN"; }

echo "── fixtures dans $WORK"
echo

# 1 — fetch d'une entrée neuve
header; row sizes.txt unverified "file://$SRC/sizes.txt" PENDING
out=$("$ROOT/scripts/fetch.sh" 2>&1); rc=$?
check "fetch d'une entrée neuve réussit" ok $rc
grep -q "PENDING" <<< "$out"
check "l'entrée neuve est signalée comme à verrouiller" ok $?
[[ -f "$RAWD/core/sizes.txt" ]]
check "le fichier est bien écrit dans raw/" ok $?

# 2 — verrouillage
"$ROOT/scripts/lock.sh" > /dev/null 2>&1
locked=$(awk -F'\t' '$1=="sizes.txt"{print $6}' "$MAN")
checkeq "le verrouillage inscrit la bonne empreinte" "$TRUE_SHA" "$locked"
status=$(awk -F'\t' '$1=="sizes.txt"{print $3}' "$MAN")
checkeq "le statut bascule en ok" "ok" "$status"
size=$(awk -F'\t' '$1=="sizes.txt"{print $7}' "$MAN")
checkeq "la taille est inscrite" "$(wc -c < "$SRC/sizes.txt" | tr -d ' ')" "$size"

# 3 — idempotence
out=$("$ROOT/scripts/fetch.sh" 2>&1); rc=$?
check "un second fetch reste valide" ok $rc
grep -qE '^\s+ok\s+sizes\.txt' <<< "$out"
check "le second fetch vérifie sans retélécharger" ok $?

# 4 — LE test : une source altérée doit être rejetée
printf 'chr1\t1\n' > "$SRC/sizes.txt"        # la source amont a changé
rm -f "$RAWD/core/sizes.txt"                  # le cache local est vidé
out=$("$ROOT/scripts/fetch.sh" 2>&1); rc=$?
check "une source altérée fait échouer le fetch" ko $rc
grep -q "MISMATCH" <<< "$out"
check "la divergence est nommée MISMATCH" ok $?
[[ ! -f "$RAWD/core/sizes.txt" ]]
check "le fichier corrompu n'est pas conservé" ok $?

# 5 — accession non résolue
header; row mystere.dat unresolved UNRESOLVED PENDING
out=$("$ROOT/scripts/fetch.sh" 2>&1); rc=$?
check "une accession non résolue fait échouer le fetch" ko $rc
grep -q "UNRESOLVED" <<< "$out"
check "l'accession manquante est nommée" ok $?

echo
printf '%d réussis, %d échoués\n' "$pass" "$fail"
((fail == 0))
