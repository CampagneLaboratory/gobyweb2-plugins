#!/bin/bash
installation_path=${RESOURCES_ARTIFACTS_ENSEMBL_API_INSTALL_DIR}
PERL5LIB=${PERL5LIB}:${installation_path}/src/bioperl-1.2.3
PERL5LIB=${PERL5LIB}:${installation_path}/src/ensembl/modules
PERL5LIB=${PERL5LIB}:${installation_path}/src/ensembl-compara/modules
PERL5LIB=${PERL5LIB}:${installation_path}/src/ensembl-variation/modules
PERL5LIB=${PERL5LIB}:${installation_path}/src/ensembl-functgenomics/modules
export PERL5LIB