
# Output directory structure

Below is the structure of the results directory. Pathways will be populated with parameters set in the [params file](../params.yml). An example of a parameters file is found in the [demo params](../demo/params__large_demo.yml).

```bash
sc_nf_diffexpression
├── differential_expression
│   ├── covariate_to_test_1
│   │   ├── [files: merged differential results]
│   │   ├── [plots: results]
│   │   ├── cell_label=cluster_1
│   │   │   ├── method=package::resolution::method___formula=cov_1pluscov_2pluscov_3
│   │   │   │   ├── [files: differential results]
│   │   │   │   ├── gsea_method-gsea_config_1=value::gsea_config_2=value::gsea_config_3=value::...::db=databases::signed
│   │   │   │   ├── gsea_method-gsea_config_1=value::gsea_config_2=value::gsea_config_3=value::...::db=databases::unsigned
│   │   │   │   ... etc. ...
│   │   │   ├── method=package::resolution::method___formula=cov_1pluscov_2pluscov_3__proportion_covs-cluster
│   │   │   ... etc. ...
│   │   ├── cell_label=cluster_2
│   │   ... etc. ...
│   ├── covariate_to_test_2
│   ... etc. ...
```
