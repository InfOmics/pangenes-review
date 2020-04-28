# pangenes-review
Pipelines used for the review of gene-oriented pan-genomic tools 

<hr>

You may want to resolve some requirements before run the pipelines
```
pip2 install ete2
pip2 install biopyton
sudo pip2 install phylo-treetime
sudo apt-get install mcl
sudo apt-get install -y raxml
sudo apt-get install -y fasttree
sudo apt-get install -y mafft
sudo apt-get install -y bioperl
sudo apt-get install -y libstatistics-linefit-perl
#curl -O https://repo.anaconda.com/archive/Anaconda3-2019.03-Linux-x86_64.sh
#bash Anaconda3-2019.03-Linux-x86_64.sh
#source /home/vbonnici/.bashrc
#conda install -c bioconda perl-statistics-distributions
sudo apt-get install -y libstatistics-distributions-perl
sudo apt-get install roary
sudo apt-get install bedtools cd-hit ncbi-blast+ mcl parallel cpanminus prank mafft fasttree
sudo cpanm -f Bio::Roary
sudo apt-get install -y libany-moose-perl
sudo apt-get install libmoose-perl
sudo apt-get install -y liblog-log4perl-perl
cpan
sudo cpanm  Array::Utils Bio::Perl Exception::Class File::Basename File::Copy File::Find::Rule File::Grep File::Path File::Slurper File::Spec File::Temp File::Which FindBin Getopt::Long Graph Graph::Writer::Dot List::Util Log::Log4perl Moose Moose::Role Text::CSV PerlIO::utf8_strict Devel::OverloadInfo Digest::MD5::File
sudo cpan File::Which
sudo cpan Array::Utils
sudo cpan File:Grep
sudo cpan Text::CSV
sudo cpan Digest::MD5::File
sudo apt-get install bedtools cd-hit blast mcl GNUparallel prank mafft fasttree

sudo cpan Parallel::ForkManager
sudo cpan Tie::Log4perl

#sudo su - -c "R -e \"install.packages('devtools')\""
#sudo su - -c "R -e \"devtools::install_github('larssnip/micropan')\""

#sudo apt purge r-base* r-recommended r-cran-*
#sudo apt autoremove
#sudo apt update
#sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
#sudo apt update
#gpg --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
#gpg -a --export E298A3A825C0D65DFD57CBB651716619E084DAB9 | sudo apt-key add -
#sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
#sudo apt update
#sudo add-apt-repository ppa:marutter/c2d4u3.5
#sudo apt-get update
#sudo apt install r-base r-base-core r-recommended
#sudo apt install r-cran-rgl r-cran-rjags r-cran-snow r-cran-ggplot2 r-cran-igraph r-cran-lme4 r-cran-rjava r-cran-devtools r-cran-roxygen2 r-cran-rjava r-cran-xlsx

sudo su - -c "R -e \"install.packages('micropan')\""
```
