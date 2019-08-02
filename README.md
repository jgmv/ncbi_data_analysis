# ncbi_data_analysis
Miscellaneous python scripts to deal with NCBI records and sequences

* **fetch_gb.py**\*: Retrieve GenBank records from NCBI from GI list.
* **get_fasta_from_gb.py**: Retrieve fasta sequences from GenBank records.
* **get_metadata_from_gb.py**: Retrieve metadata from GenBank records.
* **search_ncbi_by_term.py**\*: Retrieve GI list from NCBI from Entrez terms.
* **fetch_PubChem_compound.py**: Retrieve PubChem records from CID list.

\* Tell NCBI who you are by replacing 'your@email' with your e-mail address. Also, respect the NCBI guidelines for posting requests (see https://www.ncbi.nlm.nih.gov/books/NBK25497): 'In order not to overload the E-utility servers, NCBI recommends that users post no more than three URL requests per second and limit large jobs to either weekends or between 9:00 PM and 5:00 AM Eastern time during weekdays. Failure to comply with this policy may result in an IP address being blocked from accessing NCBI. If NCBI blocks an IP address, service will not be restored unless the developers of the software accessing the E-utilities register values of the tool and email parameters with NCBI.' Note that Biopython tools intrinsically respect the three posts per second frequency.
