SHELL := /usr/bin/env bash
HERE  := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
FETCH := $(HERE)/scripts/fetch.sh
LOCK  := $(HERE)/scripts/lock.sh
PY    := $(if $(wildcard $(HERE)/.venv/bin/python),$(HERE)/.venv/bin/python,python3)
GENO  := PYTHONPATH=$(HERE) $(PY) -m geno_pipeline

.PHONY: help data data-core data-seq data-annot data-hic data-nuclear data-clinical data-all data-lock verify selftest status build-store query tracks bench test hic-validate recon nucleus ensemble fine fine-tp53 g3d g3d-bench fiche linked search

help:
	@echo "Geno — pipeline de données"
	@echo
	@echo "  make data            core + annot (les jeux légers, à faire en premier)"
	@echo "  make data-core       tailles, bandes, trous, chaînes de conversion"
	@echo "  make data-seq        génomes de référence en 2bit (~1,6 Go)"
	@echo "  make data-annot      GENCODE, cCREs, ChromHMM, CTCF, cohésine"
	@echo "  make data-hic        Hi-C et Micro-C (plusieurs dizaines de Go)"
	@echo "  make data-nuclear    LADs, TSA-seq, Repli-seq"
	@echo "  make data-clinical   ClinVar, GWAS, contrainte, ontologies"
	@echo "  make data-all        tout"
	@echo
	@echo "  make data-lock       inscrit les empreintes des fichiers récupérés"
	@echo "  make verify          recontrôle les empreintes déjà verrouillées"
	@echo "  make selftest        prouve la vérification sans réseau"
	@echo "  make status          état du manifeste"
	@echo
	@echo "Socle 1D"
	@echo "  make build-store     construit le magasin d'intervalles depuis les fixtures"
	@echo "  make query Q=…       interroge une région, ex. Q=chr7:5,527,000-5,530,600"
	@echo "  make tracks          liste les pistes du magasin"
	@echo "  make bench [N=…]     latence de requête à l'échelle réelle"
	@echo "  make test            suite de tests"
	@echo
	@echo "Conformation"
	@echo "  make hic-validate    plante une structure Hi-C connue et valide les callers"
	@echo "  make recon           balaie l'exposant contact vers distance"
	@echo "  make nucleus [N=…]   construit un noyau diploïde complet"
	@echo "  make ensemble [N=…]  N repliements du même génome (défaut 200, ~40 min)"
	@echo "                       reprend un magasin existant ; DAMID=… ajoute la corrélation"
	@echo "  make fine            extrusion de boucles sous le TAD, et son raccord au noyau"
	@echo "                       CTCF=motifs.bed utilise des sites mesurés au lieu de plantés"
	@echo "  make fine-tp53       la même chose autour de TP53 (chr17:5,67–9,67 Mb), ses propres graines"
	@echo
	@echo "Format et streaming"
	@echo "  make g3d             écrit data/g3d/gm12878.g3d depuis les ensembles S7 et S8,"
	@echo "                       la séquence hg38 des régions fines et leurs cibles GENCODE"
	@echo "                       (et gm12878-noyau.g3d, sans aperçu, pour comparaison)"
	@echo "                       COMPARTMENTS=e1.bedGraph affiche un E1 mesuré"
	@echo "  make g3d-bench       octets écrits de 6 000 à 600 000 billes (cache des noyaux)"
	@echo
	@echo "Fiches (semaine 13)"
	@echo "  make fiche           écrit data/g3d/gm12878.fiche — une fiche par bille, sources réelles"
	@echo "  make linked          carte de contacts du modèle (gm12878.contacts) + pistes igv — semaine 14"
	@echo "  make search          index de recherche : gènes GENCODE, rsID ClinVar et GWAS (gm12878.search) — semaine 16"
	@echo "                       du manifeste — puis la vérifie (complétude, exactitude, faits)"

data:           ; @TIERS="core annot" $(FETCH)
data-core:      ; @TIERS="core" $(FETCH)
data-seq:       ; @TIERS="seq" $(FETCH)
data-annot:     ; @TIERS="annot" $(FETCH)
data-hic:       ; @TIERS="hic" $(FETCH)
data-nuclear:   ; @TIERS="nuclear" $(FETCH)
data-clinical:  ; @TIERS="clinical" $(FETCH)
data-all:       ; @TIERS="all" $(FETCH)

data-lock:      ; @$(LOCK)
verify:         ; @TIERS="all" $(FETCH)
selftest:       ; @$(HERE)/scripts/selftest.sh

build-store:    ; @$(GENO) build
tracks:         ; @$(GENO) tracks
bench:          ; @$(GENO) bench $(if $(N),--n $(N),)
test:           ; @cd $(HERE) && PYTHONPATH=$(HERE) $(PY) -m pytest
hic-validate:   ; @$(GENO) hic $(if $(N),--bins $(N),)
recon:          ; @$(GENO) recon $(if $(N),--n $(N),)
nucleus:        ; @$(GENO) nucleus $(if $(N),--bp-per-bead $(N),)
ensemble:       ; @$(GENO) ensemble $(if $(N),--n $(N),) $(if $(DAMID),--damid $(DAMID),)
g3d:            ; @$(GENO) g3d build $(if $(COMPARTMENTS),--compartments $(COMPARTMENTS),) && $(GENO) g3d build --no-preview $(if $(COMPARTMENTS),--compartments $(COMPARTMENTS),) --out $(HERE)/data/g3d/gm12878-noyau.g3d >/dev/null && $(GENO) g3d inspect
g3d-bench:      ; @$(GENO) g3d bench
fiche:          ; @$(GENO) fiche build && $(GENO) fiche check
linked:         ; @$(GENO) linked build && $(GENO) linked check
search:         ; @$(GENO) search build && $(GENO) search check
fine:           ; @$(GENO) fine $(if $(N),--bp-per-bead $(N),) $(if $(CTCF),--ctcf $(CTCF),)
# Semaine 12. D'autres graines que la région ACTB : avec les mêmes, le 1D (numpy pur, graines
# fixes) planterait exactement les mêmes sites aux mêmes rangs, et rendrait les mêmes boucles.
fine-tp53:      ; @$(GENO) fine --chrom chr17 --start 5670000 --end 9670000 --site-seed 17 --first-seed 17000 --out $(HERE)/data/fine/tp53.npz $(if $(CTCF),--ctcf $(CTCF),)

query:
	@test -n "$(Q)" || { echo "usage : make query Q=chr7:5,527,000-5,530,600" >&2; exit 2; }
	@$(GENO) query "$(Q)" --time

status:
	@awk -F'\t' '\
	  /^#/ || $$1=="id" || NF<2 {next} \
	  {n[$$3]++; t[$$2]++; total++} \
	  END { \
	    printf "%d entrées\n\n", total; \
	    printf "  par statut\n"; for (k in n) printf "    %-12s %d\n", k, n[k]; \
	    printf "\n  par tier\n";  for (k in t) printf "    %-12s %d\n", k, t[k]; \
	  }' $(HERE)/data/manifest.tsv
