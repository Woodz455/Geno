# Délègue au pipeline. `make data` depuis la racine du dépôt fait ce qu'annonce la feuille de route.
.PHONY: help data data-core data-seq data-annot data-hic data-nuclear data-clinical data-all data-lock verify selftest status build-store query tracks bench test

help data data-core data-seq data-annot data-hic data-nuclear data-clinical data-all data-lock verify selftest status build-store query tracks bench test:
	@$(MAKE) --no-print-directory -C pipeline $@ Q="$(Q)" N="$(N)"
