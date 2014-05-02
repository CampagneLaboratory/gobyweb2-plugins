gobyweb2-plugins
================

Plugin repository for GobyWeb 2.2+. This project hosts plugins for GobyWeb versions 2.2+. 
These plugins are designed to install automatically when run on the compute nodes of the grid.
For more details about GobyWeb, see http://gobyweb.campagnelab.org (web-app) and 
http://campagnelab.org/software/gobyweb/plugins-sdk/ (plugin software development kit to use 
plugins from the command line and facilitate their development).

Please feel free to submit new plugins to this repository as pull requests. Plugins that require
specific software or data not found on a linux node by default must be implemented to provide 
auto-installation of these artifacts.

Note that since GobyWeb 2.3, you need to checkout the plugins-SDK branch of this repository. All 
new developments are performed on this branch. Changes will be merged with master after we have 
upgraded the GobyWeb CTSC production instance that depends on master.

