# Description of EDC
The Electronic Data Capture (EDC) System used in this study is called Labmatrix,
developed by BioFortis. Labmatrix is agnostic to research domain and is designed
to manage data in a secure, centralized repository. The data servers are housed
within NHGRI's secure network, behind NIH's firewall, so that only users with an
NIH account, within the NIH network or logged in through the VPN, can access it.

Labmatrix provides comprehensive integration with other applications and
database systems, e.g. electronic medical records, clinical trial management
systems, LIMS, and other research applications such as genomic or proteomic
databases. This enables flexible user-driven electronic data capture forms
without additional programming or customization. Labmatrix also provides
flexible and fine-grained user access control and requires minimal training. It
conforms to regulatory compliance standards (e.g. HIPAA, 21CFR.11) and secures
PHI. Labmatrix also stores information regarding record modification, providing
an audit trail and electronic signature features. It allows for data exporting
as commonly used formats (e.g. CSV, XLSX) to enable data analysis in any major
statistical software.

The version currently used is 6.3.3 (Revision: 23600). 

# Additional data storage
Raw data that cannot be stored in our EDC (e.g. imaging or genomics data) are
stored in the lab's network drive. Access to the drive is restricted to members
of the lab through NHGRI's Information Technology team. Similarly to Labmatrix,
the drive resides within NHGRI's secure network, behind NIH's firewall, so that
only users with an NIH account, within the NIH network or logged in through the
VPN, can access it.

Different types of data are stored in separate directories in the shared drive,
indexed by randomized identifiers. The link between those IDs and patient
information is stored in Labmatrix and in a hard-copy stored in a secure
location in the lab.

# Data management design
The data for this study is largely stored through Custom Data forms in
Labmatrix. Each form can have many fields, and Labmatrix provides several
different types of fields, which include: integers, decimals, dates, times,
notes, Yes/No, and lists. 

For questionnaires, each form mimics the question structure in the paper copy to
facilitate data entry. Some questionnaires might have aggregated scores in the
end, which is implemented in Labmatrix through formulas written in Javascript
code. When the researcher finishes manually entering the data, the aggregated
scores are calculated and saved with the corresponding record.

Neuropsychiatric tests are stored in a similar fashion to questionnaires, where
each form stores the results of a specific test. Manual test scoring is needed,
which is dependent on each test, and those results are then entered into
Labmatrix.

Only metadata about neuroimaging scans is saved to Labmatrix. Those include the
scan identifier, type of scan, date, scanner, and any notes. The actual neural
data is stored in the lab network drive and indexed by the identifier. Results
of visual quality control for each scan are also stored in Labmatrix.

Biosamples data require a more complex form of storage in Labmatrix, but it is
still implemented through a Custom Data form. The lab tracks each movement from
the date when a sample is collected until genomic data is returned to us about
the sample. That requires storing not only biosample type (e.g. Blood, Saliva,
Blood DNA, SNP data, etc), but also the sample status (e.g. in inventory,
shipped for DNA extraction, etc), physical location, and day in / day out (e.g.
date of collection, and date of shipment for extraction). Similarly to
neuroimaging, the actual genomics data are stored in the lab shared network
drive, indexed by their unique IDs.

The only electronic data file stored in Labmatrix is the consent form (scanned
PDF). The Custom Data form to store the file also includes fields to store
participants' answers about data sharing that are part of the consent process.

Finally, the data collected during clinical interviews, as well as data acquired
by computerized instruments (e.g. CPT, response tp behavioral tasks while in the
MRI scanner) are backed up regularly and stored in the shared drive.

# Data integrity checks
All fields in Labmatrix Custom Data forms share the option of whether they are
required to be filled in or not. Then, Labmatrix provides a Javascript interface
that enables the input of code to check any entries. For example, valid-range
checks, missing values, and formulas for automatic filling in fields based on
other values can be set up this way.

Labmatrix also periodically runs data integrity checks with respect to subject
demographics. Those tests validate the demographics data entered in Labmatrix
against the data entered in the Clinical Center database during the admissions
process. Incongruencies are highlighted to be later manually fixed by
researchers in the lab.

Data are also exported from Labmatrix on a monthly basis to run offline scripts
in R and Python that alert for multivariate and cross-module errors and ensure
data consistency. Examples of such errors can be the same biosample ID assigned
to different participants (which is later corrected based on hardcopies of ID
assignments), date mismatch between visits and data collection, and sample
history monitoring error (for example, a biosample DNA extraction date cannot
happen before the sample collection date).

# Data recovery
Labmatrix keeps an audit trail and regular backups in case the data need to be
restored to a specific time point. The background database is housed in
redundant servers under the NIH firewall.

The NHGRI's IT team performs regular backups of the shared drive to mitigate
data loss.


# Data priviledges and roles

Any person in the lab, prior to being granted permission to use Labmatrix, must
be trained on how to use the EDC by Labmatrix staff. After this general
Labmatrix training, each new person also goes through a more specific Labmatrix
training with the data manager to ensure understanding of the specific design
for data storage in the lab.

All lab members have permissions to enter and query data in Labmatrix after
initial training, but only the data manager can remove records.

For the most part, the data are based on source documentation maintained at the
clinical site, either in digital or hard copy. For cases where the EDC system
stores only metadata (scans and biosamples, for example), data are stored and
organized in the network drive by the data manager.


# Data queries and export
All data queries in Labmatrix are generated through a tool called Qiagram. It is
a visual interface to construct a SQL query on the background, by creating set
operations (intersection, union, filtering, etc) transparent to the user. Data
queries can be performed by all members of the lab, as they do not alter the
data in Labmatrix. Results are then exported from Labmatrix for any sort of data
analysis. It supports CSV, tab-delimited, XML, and Excel formats, which are
easily imported by any major data analysis tool.

Sometimes data are also exported to aid editing, as there are operations that
are easier to be performed in a matrix format (i.e. in Excel), and then
re-uploaded to Labmatrix. Labmatrix has the cabalitity of "checking out" records
for editing, so that when re-uploaded only the appropriate records are changed.
