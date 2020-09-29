# fhir_converter
The purpose of this program to convert from patients varients into FHIR json to the format that required by Pharmcat  
This file was written by python, and the version is Python 3.7.6

In order to use this fhir_converter.py, all files need to be in a same directory. There are two files that required for running this program.
- PharmCAT conversion table
  Contains the conversion between b37 SPDI to b38 SPDI
- PharmCat Template table

Also, at the end of this program, it requires PharmCAT installed to achieve a final report from PharmCAT. 
  The installation and more informations about PharmCAT : https://github.com/PharmGKB/PharmCAT/wiki

To run the program, use command line on a terminal: python3 fhir_converter.py FHIR_json_file_name

As an example, if the FHIR json file name is "NA162.CYP2C19.fhir.json", the command line should be:
```
$ python3 fhir_converter.py NA162.CYP2C19.fhir.json
```

After running this command line, new file : NA162.CYP2C19.fhir.json.vcf will be created in the same directory, adding the exact vcf file name after the following command line
java -jar pharmcat-0.7.0-all.jar -vcf and -o output which will create new directory called output. 

```
$ java -jar pharmcat-0.7.0-all.jar -vcf NA162.CYP2C19.fhir.json.vcf -o output
```

If all the process are successfully completed, PharmCAT report in HTML format will be in the output directory under same name as the previous vcf file. 

Ex: NA162.CYP2C19.fhir.json.vcf.html
