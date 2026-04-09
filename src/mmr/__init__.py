"""MMR deficiency classifier modules.

Original author: Ayan Choudhury
Integrated into the project pipeline by Keith Coughlin (paths, config,
import conventions adapted to match the rest of src/).

Trains a within-COAD classifier that predicts mismatch repair (MMR) protein
loss status (dMMR vs pMMR) from gene expression. Companion to the omics-only
biomarker scoring pipeline; addresses a within-cancer subtype question
relevant to immune checkpoint inhibitor response.
"""
