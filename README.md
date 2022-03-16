# Library of RNA Bioinformatics Algorithms
This library provides basic algorithms in RNA bioinformatics.
Available current algorithms are the [gamma-centroid structure estimator](https://academic.oup.com/bioinformatics/article/25/4/465/248190), [McCaskill's algorithm](https://onlinelibrary.wiley.com/doi/abs/10.1002/bip.360290621), and [Durbin's algorithm](http://www.mcb111.org/w06/durbin_book.pdf).
To model RNA secondary structure scoring, [the Turner 2004 model](https://rna.urmc.rochester.edu/NNDB/) and [the CONTRAfold v2.02 model](http://contra.stanford.edu/contrafold/) are available.
I adopt [the CONTRAlign v2.01 model](http://contra.stanford.edu/contralign) to model RNA sequence alignment scoring.
More descriptions of these two scoring models are [here](https://github.com/heartsh/rna-ss-params).
The gamma-centroid structure estimator of this repository is compared to that of [CentroidFold](https://github.com/satoken/centroid-rna-package) regarding secondary structure prediction accuracy:

|![Gamma-centroid estimator performance of different models and different implementations](./assets/images_fixed/fig_1.png)
|:-:|

# Author
[Heartsh](https://github.com/heartsh)

# License
Copyright (c) 2018 Heartsh  
Licensed under [the MIT license](http://opensource.org/licenses/MIT).
