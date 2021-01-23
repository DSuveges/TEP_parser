# TEP retrieval and parser

This script retrieves [TEP]() (target enabling packages) data from the [SGC]() website enrich with target annotation using Ensembl REST api then saves the data in a gzipped json file.

## Usage

```bash
python TEP_retrieve.py --output <outputFile> --logFile <logfile name>
```

**Where**:

* `outputFile`: gzipped JSON object following the structure shown below
* `logfile`: file intowhich the logs are saved. Optional. If not provided logs are written to the standard error

## Output

The parsed data saved in this format accomodating the needs of the OT frontend:

```json
{
  "ENSG00000143379": {
    "id": "ENSG00000143379",
    "symbol": "SETDB1",
    "link": "https://www.thesgc.org/tep/SETDB1",
    "disease": "Oncology",
    "uniprot_id": "Q15047"
  },
  "ENSG00000167258": {
    "id": "ENSG00000167258",
    "symbol": "CDK12",
    "link": "https://www.thesgc.org/tep/CDK12",
    "disease": "Oncology",
    "uniprot_id": "Q9NYV4"
  }
}
```



