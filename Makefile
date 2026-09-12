# Délègue au pipeline. `make data` depuis la racine du dépôt fait ce qu'annonce la feuille de route.
.PHONY: help data data-core data-seq data-annot data-hic data-nuclear data-clinical data-all data-lock verify selftest status

help data data-core data-seq data-annot data-hic data-nuclear data-clinical data-all data-lock verify selftest status:
	@$(MAKE) --no-print-directory -C pipeline $@
