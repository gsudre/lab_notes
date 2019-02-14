# 2019-02-14 09:12:45

Let's keep track of how we're running the FDT pipeline on our own data. First,
let's do it only for the people that have genotype data, as they'l likely be
used for gene expression project.

First, copy everything to the cluster so we can batch convert the data. We can
create a Globus script for that, getting the DICOM folders straight from the
README file to avoid any interactions with the previous DTI preprocessing.

