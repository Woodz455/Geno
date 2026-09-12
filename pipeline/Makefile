SHELL := /usr/bin/env bash
HERE  := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
FETCH := $(HERE)/scripts/fetch.sh
LOCK  := $(HERE)/scripts/lock.sh

.PHONY: help data data-core data-seq data-annot data-hic data-nuclear data-clinical data-all data-lock verify selftest status

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

status:
	@awk -F'\t' '\
	  /^#/ || $$1=="id" || NF<2 {next} \
	  {n[$$3]++; t[$$2]++; total++} \
	  END { \
	    printf "%d entrées\n\n", total; \
	    printf "  par statut\n"; for (k in n) printf "    %-12s %d\n", k, n[k]; \
	    printf "\n  par tier\n";  for (k in t) printf "    %-12s %d\n", k, t[k]; \
	  }' $(HERE)/data/manifest.tsv
