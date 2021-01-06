---
title: "Automatic and Interpretable Model for Periodontitis Diagnosis in Panoramic Radiographs"
date: 2020-03-14
authors: [
    "Haoyang Li", 
    "Juexiao Zhou",
    "**Yi Zhou**",
    "Jieyu Chen",
    "Feng Gao",
    "Ying Xu",
    "Xin Gao"
    ]  # Markdown supported
featured: false
journal: "MICCAI"
doi: 
pdf_url: 
cite_url: 
---

Periodontitis is a prevalent and irreversible chronic inflammatory disease both in developed and developing countries, and affects about 20% - 50% of the global population. The tool for automatically diagnosing periodontitis is highly demanded to screen at-risk people for periodontitis and its early detection could prevent the onset of tooth loss, especially in local community and health care settings with limited dental professionals. In the medical field, doctors need to understand and trust the decisions made by computational models and proposing interpretable machine learning models is crucial for disease diagnosis. Based on these considerations, we propose an interpretable machine learning method called Deetal-Perio to predict the severity degree of periodontitis in dental panoramic radiographs. In our method, alveolar bone loss (ABL), the clinical hallmark for periodontitis diagnosis, could be interpreted as the key feature. To calculate ABL, we also propose a method for teeth numbering and segmentation. First, Deetal-Perio segments and indexes the individual tooth via Mask R-CNN combined with a novel calibration method. Next, Deetal-Perio segments the contour of the alveolar bone and calculates a ratio for each individual tooth to represent ABL. Finally, Deetal-Perio predicts the severity degree of periodontitis given the ratios of all the teeth. The entire architecture could not only outperform state-of-the-art methods and show robustness on two data sets in both periodontitis prediction, and teeth numbering and segmentation tasks, but also be interpretable for doctors to understand the reason why Deetal-Perio works so well.
