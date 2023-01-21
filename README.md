# Phage-Host Interaction Prediction

## Description
**ABSTRACT**: With the growing interest in using phages to combat antimicrobial resistance, computational methods for predicting phage-host interactions have been explored to help shortlist candidate phages. While existing systems have been successful in integrating multiple features to improve performance, most consider entire proteomes and rely on manual feature engineering, which poses difficulty in selecting the most informative sequence properties to serve as input to the model. In this paper, we focused on the phages' receptor-binding proteins, which are known to be the key machinery for host recognition, and explored different protein language models to automatically encode these protein sequences into meaningful dense embeddings without the need for additional alignment or structural information. Our experiments showed that the use of embeddings of receptor-binding proteins presents improvements over handcrafted genomic and protein sequence features for phage-host interaction prediction. The highest performance was obtained using the transformer-based protein language model ProtT5, resulting in a 3% to 4% increase in the weighted F1 scores across different thresholds for prediction confidence.

## Project Structure

## Authors
- **Mark Edward M. Gonzales** <br>
  mark_gonzales@dlsu.edu.ph 
 
- **Ms. Jennifer C. Ureta** <br>
  jennifer.ureta@dlsu.edu.ph 
  
- **Dr. Anish M.S. Shrestha** <br>
  anish.shrestha@dlsu.edu.ph

This is a research project under the [Bioinformatics Laboratory](https://bioinfodlsu.com/), Advanced Research Institute for Informatics, Computing and Networking, De La Salle University, Philippines.
  
