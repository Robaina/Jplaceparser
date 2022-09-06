![logo](assets/logo.png)

Python tools to parse and manipulate [JPlace files](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0031009), the format of [Evolutionary Placement](https://arxiv.org/abs/0911.2852) results.

[![DOI](https://zenodo.org/badge/428328886.svg)](https://zenodo.org/badge/latestdoi/428328886)

<br>

JplaceParser allows not only reading Jplace files as python dictionaries but also manipulating field values and exporting back to jplace files. It contains a number of filters to remove placements that do not satisfy quality criteria. Currently, placements can be filtered by three criteria involving the [pendant and distal](https://github.com/lczech/gappa/wiki/Subcommand:-assign#automatic-ratio-example) placement length as well as the phylogenetic tree diameter:

* Filter by maximum pendant length
* Filter by maximum pendant to distal length ratio
* Filter by maximum pendant to tree diameter ratio

This is an ongoing project!

## Installation
1. ```pip install jplaceparser```

or

2. Git clone project to local directory.

   In terminal navigate to directory and enter: ```python setup.py install```


```python
from jplaceparser import JplaceParser


jplace = JplaceParser.fromJplaceFile("examples/example.jplace")
jplace
```





<table>
    <tr>
        <td><strong>Number of Placements</strong></td><td>19</td>
    </tr><tr>
        <td><strong>Fields</strong></td><td>edge_num, likelihood, like_weight_ratio, distal_length, pendant_length</td>
    </tr><tr>
        <td><strong>JplaceParser version</strong></td><td>0.0.1</td>
    </tr><tr>
        <td><strong>Author</strong></td><td>Semidán Robaina Estévez, 2022</td>
    </tr>
</table>





```python
filtered_jplace = jplace.filterByMaxPendantToTreeDiameterRatio(
    max_pendant_ratio=0.001
)

filtered_jplace.writeToFile("examples/filtered_example.jplace")

filtered_jplace
```

    Filtering placements for tree diameter: 4.519636416






<table>
    <tr>
        <td><strong>Number of Placements</strong></td><td>9</td>
    </tr><tr>
        <td><strong>Fields</strong></td><td>edge_num, likelihood, like_weight_ratio, distal_length, pendant_length</td>
    </tr><tr>
        <td><strong>JplaceParser version</strong></td><td>0.0.1</td>
    </tr><tr>
        <td><strong>Author</strong></td><td>Semidán Robaina Estévez, 2022</td>
    </tr>
</table>



