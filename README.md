TODO: Create kubeflow pipeline for the entire process.
TODO: Finalize README.md

Expert Knowledge:

1) Bases can be A, C, G, T, N (https://pubs.acs.org/doi/abs/10.1021/ar200257x)
2) For each chromosome, the length in the metadata should be larger than the value in column POS for that chromosome. (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.39/)

Statistical Analysis:

Check the google doc file and also the underscored text in the VCFv4.2 pdf file.

Adversarial Evaluation:

Split the data into original and synthetic samples.
Let's set an 80% overlap between the samples in these two sets.

Get the output csv for the original samples and then for the synthetic samples.

Then treat each column as a different label (a different task). Then each task is a multi-class classification problem and we can calculate metrics such as accuracy, precision, recall, F1 score, etc.

## 📜 License & Usage

All rights reserved by MetaMinds Innovations.
