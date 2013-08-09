GobyWeb2 Plugins Change log
===========================

## Towards release 2.2.2 (future)
1. ENSEMBL_ANNOTATIONS resource:  Biomart.groovy script now restarts dropped connections up to three times. Advance version number to 2.1.2.

## Release 2.2.1 (Aug 9 2013)

1. Fix VEP plugin for new version of Variant Effect Predictor (v72). Upgrade links for ensembl API distribution. Fetch and decompress ensembl-tools archive because VEP is now packaged in this file (would have been nice to see this listed as a change in the VEP/Ensembl API change log.
2. ENSEMBL_ANNOTATIONS resource: Biomart.groovy script. Add filter by chromosome for new annotation kinds. This attribute forces download chromosome by chromosome, which is necessary with Biomart for large tables.
3. ENSEMBL_ANNOTATIONS resource: Add pre-packed cache to speed up biomart annotation installs. The cached file is written and read from ~gobyweb/url-cache. Edit the plugin as needed to suit your particular local installation. 
4. upgrade GOBYWEB_SERVER_SIDE/global_goby.jar to Goby 2.3.2 pre-release
5. Expose BWA aln options in BWA_GOBY_ARTIFACT plugin.
6. PLAST (parallel last) plugin: Extend PLAST to support paired end reads. New version of Goby that supports --paired argument in run-parallel. 
