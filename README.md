# Density-Classification-and-Microclacifition-Detection-in-Mammogram-Diagnosis

Implemented methods to automate the classification of breast densities and detection of microcalcification regions

Preprocessing pipeline included left-right flipping, application of a 3x3 median filter, binary mask generation, uniform image resizing, and removal of pectoral muscle
  
  
Analyzed texture via Gabor filtering and extracted 23 features
  
  
Utilized a morphological based approach for microcalcification segmentation
  
  
Classified density using quadratic discriminant analysis with
three-fold random cross-validation (90-10 split)
