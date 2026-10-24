"""Recherche (semaine 16) : un gène, une position, un rsID → un intervalle hg38.

`query` lit la requête telle qu'on la tape ; `build` écrit l'index (`gm12878.search`) qu'elle
interroge — les gènes de GENCODE et les rsID de ClinVar et du catalogue GWAS, rien d'autre ;
`check` le recompte depuis les fichiers bruts. Le lecteur TypeScript
(`packages/g3d/src/search.ts`) applique la même grammaire au même fichier.
"""
