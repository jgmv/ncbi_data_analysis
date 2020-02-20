# ncbi_data_analysis
Miscellaneous python scripts to deal with NCBI records and sequences

### draw_chem.py
Draws chemical structures from a file containing CID and SMILES values (in each row, space-separated).

### fetch_gb.py[^1]
Retrieves GenBank records from NCBI from GI list.

### fetch_PubChem_compound.py
Retrieves PubChem records from CID list.

### get_fasta_from_gb.py
Retrieves fasta sequences from GenBank records.

### get_metadata_from_BioSample.py
Retrieves metadata from BioSample records' summary as a table.

### get_metadata_from_gb.py
Retrieves metadata from GenBank records.

### get_prot_from_gb.pyRetrieves protein sequences metadata from GenBank records.

### search_ncbi_by_term.py[^1]
Retrieve GI list from NCBI from Entrez terms.

[^1] Tell NCBI who you are by stating your e-mail address using `-email <your@email>`. Also, respect the NCBI guidelines for posting requests (see [https://www.ncbi.nlm.nih.gov/books/NBK25497](https://www.ncbi.nlm.nih.gov/books/NBK25497)): 'In order not to overload the E-utility servers, NCBI recommends that users post no more than three URL requests per second and limit large jobs to either weekends or between 9:00 PM and 5:00 AM Eastern time during weekdays. Failure to comply with this policy  may result in an IP address being blocked from accessing NCBI. If NCBI blocks an IP address, service will not be restored unless the developers of the software accessing the E-utilities register values of the tool and email parameters with NCBI.' Note that Biopython tools intrinsically respect the three posts per second frequency.
