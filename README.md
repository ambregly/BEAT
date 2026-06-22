# BEAT — Comparaison des fusions kmer vs vizome (BEAT AML)

Cette procédure compare les fusions de gènes détectées par une approche **kmer** aux
fusions de référence de **vizome / BEAT AML**, calcule les vrais positifs (TP),
faux positifs (FP) et faux négatifs (FN), la **précision** et le **recall**, et produit
des tableaux et des graphiques.

---

## 1. Installation

Python 3 suffit. Installer les dépendances une fois :

```bash
pip install -r requirements.txt
```

(matplotlib, matplotlib-venn, pandas, openpyxl, scipy, numpy)

Un script R équivalent pour les visualisations est aussi fourni (`visualize_fusions.R`,
dépendances installées automatiquement au premier lancement).

---

## 2. Format des fichiers d'entrée

### `BEAT_AML.csv` (la référence « vizome »)
CSV avec en-tête. Colonnes utilisées :
`SampleID, left_gene, left_chr, right_gene, right_chr, fusion_index`
(+ `junction_read_count`, `spanning_frag_count`, `ffpm` pour les comptages).
Une ligne = **une fusion** (identifiée par `SampleID` + `fusion_index`).
Exemple de `SampleID` : `BA2409`.

### `kmer.tsv` (les détections kmer)
TSV **sans en-tête**, 3 colonnes :

| Colonne | Contenu | Exemple |
|---|---|---|
| 1 | `left_gene_left_chr_pos_pos_right_gene_right_chr_pos_pos_idx\|idx\|idx` | `ELOVL5_chr6_53213614_53213640_PTP4A1_chr6_64286340_64286365_2744\|3218\|3604` |
| 2 | `SampleID` avec un suffixe `R` | `BA2409R` |
| 3 | comptage kmer | `7.0` |

Les `idx` de la colonne 1 sont les `fusion_index` de BEAT AML associés à cette signature.
Le suffixe `R` du `SampleID` est retiré automatiquement (`BA2409R` → `BA2409`).
Le chromosome est harmonisé (`chr6` ↔ `6`).

---

## 3. Comment une fusion est comptée (logique de matching)

Une fusion est identifiée par
`(SampleID, left_gene, left_chr, right_gene, right_chr, fusion_index)`.

Mode par défaut **`row`** (recommandé) : l'unité comptée est **la ligne kmer**.
- Ses `fusion_index` cherchent leur correspondance (échantillon + index) dans BEAT AML.
- Dès qu'**un** index correspond → **1 vrai positif** (les autres index de la ligne sont ignorés).
- Si **aucun** ne correspond → **1 seul faux positif**.
- **FN** = fusions de BEAT AML jamais retrouvées.

Autres modes (`--match`) :
- `index` : chaque couple (échantillon, index) est compté individuellement.
- `genepair` : correspondance par (échantillon + paire de gènes), `fusion_index` ignoré.

`precision = TP/(TP+FP)`, `recall = TP/(TP+FN)`, `F1 = 2·P·R/(P+R)`.

---

## 4. Usage rapide : tout en une commande

Le script `run_pipeline.py` enchaîne les 3 étapes (comparaison → tableaux/Venn → graphes).

```bash
python3 run_pipeline.py --beat BEAT_AML2.csv --kmer kmer2.tsv \
    --label mon_run --common-samples --match row --outdir analyses
```

Sorties créées sous `analyses/mon_run/` :

```
resultats/
    vrais_positifs.tsv      faux_positifs.tsv      faux_negatifs.tsv
figures/
    confusion_par_fusion.csv / .pdf            # 1 ligne par fusion
    confusion_par_echantillon.csv / .pdf       # 1 ligne par échantillon
    confusion.xlsx                             # 2 onglets (fusion / échantillon)
    venn_kmer_vizome.png                       # Venn kmer vs vizome
    scatter_kmer_vs_vizome.png                 # comptage kmer vs vizome (TP)
    hist_f1_par_fusion.png                     # répartition des F1
    faux_negatifs_vizome.csv                   # liste FN enrichie du comptage
    tp_comptages.csv                           # données du nuage de points
```

### Plusieurs jeux / paramètres d'un coup (manifeste)

Créer un fichier TSV (voir `jeux_exemple.tsv`) avec une ligne par run :

```
label	beat	kmer	match	common_samples	rows_per_page	vizome_count
run_all	BEAT_AML2.csv	kmer2.tsv	row	no	30	junction_read_count
run_common	BEAT_AML2.csv	kmer2.tsv	row	yes	30	junction_read_count
run_index	BEAT_AML2.csv	kmer2.tsv	index	yes	30	junction_read_count
```

```bash
python3 run_pipeline.py --manifest jeux_exemple.tsv --outdir analyses
```

Chaque run a son propre dossier `analyses/<label>/`, ce qui permet de comparer
facilement les paramètres côte à côte.

---

## 5. Les paramètres avec lesquels jouer

| Paramètre | Valeurs | Effet |
|---|---|---|
| `--match` | `row` (défaut), `index`, `genepair` | unité de comptage des TP/FP (voir §3) |
| `--common-samples` | présent/absent | ne garde que les échantillons présents **dans les deux** fichiers (réduit les FN inutiles dus aux échantillons non partagés ; améliore le recall) |
| `--vizome-count` | `junction_read_count` (défaut), `spanning_frag_count`, `junction_plus_spanning`, `ffpm` | colonne BEAT AML servant d'axe Y du nuage de points |
| `--rows-per-page` | entier (défaut 40) | nombre de lignes par page dans les tableaux PDF |
| `--outdir` | chemin | dossier racine des sorties |
| `--label` | texte | nom du sous-dossier du run |

---

## 6. Lancer les étapes séparément (optionnel)

```bash
# a) comparaison -> TP/FP/FN
python3 compare_fusions.py BEAT_AML2.csv kmer2.tsv --outdir resultats \
    --match row --common-samples

# b) tableaux + Venn (à partir du dossier resultats)
python3 visualize_fusions.py resultats --outdir figures --rows-per-page 30
#   ou en R :
Rscript visualize_fusions.R resultats figures 30

# c) graphes supplémentaires (scatter, F1, liste FN) à partir des fichiers sources
python3 extra_plots.py BEAT_AML2.csv kmer2.tsv --outdir figures \
    --match row --common-samples --vizome-count junction_read_count
```

---

## 7. Fichiers du dépôt

| Fichier | Rôle |
|---|---|
| `compare_fusions.py` | comparaison kmer vs BEAT AML → TP/FP/FN, précision, recall |
| `visualize_fusions.py` | tableaux de confusion (CSV/Excel/PDF) + diagramme de Venn |
| `visualize_fusions.R` | équivalent R des visualisations |
| `extra_plots.py` | nuage de points kmer vs vizome, histogramme F1, liste FN enrichie |
| `run_pipeline.py` | orchestrateur (un jeu ou plusieurs via manifeste) |
| `jeux_exemple.tsv` | manifeste modèle |
| `requirements.txt` | dépendances Python |

---

## 8. Note sur « combien de fusions »

- **Nombre de fusions** = nombre de lignes (couples échantillon + fusion_index).
- **Nombre d'échantillons** = `SampleID` distincts.

Un échantillon porte plusieurs fusions (souvent ~4). Donc, par exemple,
« 356 échantillons » et « 1501 fusions » ne sont pas contradictoires.
Dans le diagramme de Venn, le cercle « kmer » = TP+FP (lignes kmer) et
le cercle « vizome » = TP+FN (lignes BEAT AML).
