program then stops and prints out a message that looks like the one below:
```
. Site = 267
. invar = -1
. scale_left = 0 scale_rght = 0
. Lk = 0 LOG(Lk) = -inf < -1.79769E+308
. rr=0.136954 p=0.250000
. rr=0.476752 p=0.250000
. rr=1.000000 p=0.250000
. rr=2.386294 p=0.250000
. pinv = 0.2
. Err in file lk.c at line 739
```
If this occurs with your data set, please use the patch biglk.patch available from the Downloads page on this site. To apply this patch, just go in the src/ directory and type 'patch -i biglk.patch'. You will then need to compile the program again. Use this patch with version 20120412 of PhyML only.

  * _19 September 2011._ Versions prior to phyml-20110919.tar.gz have an important bug affecting PhyTime. **Priors on node ages are not recorded properly** (except for the root node). Analysis conducted with versions of PhyML older that this one need to be re-run. Sorry for the inconvenience.

  * _04 March 2011._ Due to a bug introduced in [revision 568](https://code.google.com/p/phyml/source/detail?r=568) (commited the 3rd of Nov. 2010), the bootstrap analysis would not complete. It has been fixed in [revision 596](https://code.google.com/p/phyml/source/detail?r=596).

  * _03 Feb 2010._ Bug found in **bootstrap calculation**. The comparison of bipartition of equal size could fail when the first character of taxon names were identical. I expect the impact of this bug to be minor.

  * _03 Feb 2010._ The calculation of **unconstrained likelihood** values ignore gap and indels. It is not a bug but a more sophisticated treatment of those characters would be more appropriate.

  * _14 Sep 2010._ The calculation of **bootstrap values with the MPI version of PhyML development version** is flawed. This bug is expected to have been introduced in the program early in 2010. It does not affect the non-development version of PhyML nor the binaries available on http://www.atgc-montpellier.fr.