# MFLAB README  (Heemstede, Netherlands, 8 Dec 2015)

`mfLab` ("Modflow Laboratory") is a n open-source scripting envioronment in Matlab
to set up, run and analyze groundwater flow and transport models of
the suite MODFLOW (different versions), MODPATH, MT3DMS and SEAWAT (all open source).

`mfLab` runs on Macs and Windows and was developed on a MacBook.pro.

`mfLab` requires `Matlab` (the Mathworks(c)), which is not open source.

`mfLab` has been available as open-source groundwater modeling environment since 2008 with
1088 commits (updates) during its development between 2008 and 2015. It waas previously available
on `http://code.google.com/p/mflab` and has been successfully used in many groundwater
modelling projects; all my previous MSc and PhD students at TUDelft have used it for their grounwater
models in their research for their theses. Also my colleagues in Waternet (www.waternet.nl)
as well as a number of people worldwide, to solve small and big groundwater problems.

It has also often been used to construct smaller models, cross sections both flat and axially symmetric, meant for
demonstration purposes as well as in courses. Many example simulation movies can be found
on YouTube (search for "groundwater" and "olsthoorn" of "mflab").

Most intensively `mflab` has been used to solve groundwater problems affected by
density flow caused by varying salinities as well as by temperature, which also includes the
effects of viscosity. Hence, most of the examples can be found in the SEAWAT folder. Further, notice
that some examples are provided that simulate fresh and saline groundwater flow using the salt water intrusion package 
SWI, which is included in MF2005.

Because Google stopped supporting their `code.google.com` site in 2015, the project had to
be moved to `github`. As there were many problems with this transfer, I had to establish a
clean version of the most recent files and create a new `mfLab` repository on Github.
Therefore, the current Github repository does not contain the the subversion history that
has been available from the early days of the development in te googlecode reposiotory, hence, a new history will be built up from
December 8 2015 on Github.

The repository also contains documentation and tutorials.

Theo N. Olsthoorn (prof.dr.ir., emer. University of Technology Delft and emer. hydrologist at Waternet)
tolsthoorn@gmail.com
Heemstede, the Netherlands

# Update from mflab (using Matlab) to mf6lab (using Python and flopy with Modflow 6)

Since 2017 I moved from using Matlab to using Python because Python is open-source and therefore students and users are no longer depend on pricy Matlab being available to them. Onn top of this, the developers of Modflow at USGS also transferred to Python and scripting to build the input of Modflow models, which which they developed `flopy`. `flopy`. For this reason I too moved to using flopy for my modeling.

I started a new follow-up project named **mf6lab** in pythhon to more easily construct groundwater models in conjunction with flopy and underlying new version of the USGS groundwater modellering code, called **MODFLOW6** or **MF6**. The repository **mf6lab** can be found on my github page.

Theo N. Olsthoorn, 2024-02-02
