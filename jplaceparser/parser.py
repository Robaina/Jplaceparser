from __future__ import annotations
from importlib import metadata
import re 
import json 
from pathlib import Path
from io import StringIO

from Bio import Phylo

meta = metadata.metadata("jplaceparser")
__version__ = meta["Version"]
__author__ = meta["Author"]


class JplaceParser():
    """
    Methods to parse jplace files, as specified in 
    https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0031009
    """
    def __init__(self, jplace_object: dict) -> None:
        self._jplace_obj = jplace_object
        self._tree_obj = next(
            Phylo.parse(
                StringIO(self.newickfyTree(self._jplace_obj['tree'])), 'newick'
                )
            )

    def _repr_html_(self):
        """This method is executed automatically by Jupyter to print html!"""
        return """
        <table>
            <tr>
                <td><strong>Number of Placements</strong></td><td>{n_plac}</td>
            </tr><tr>
                <td><strong>Fields</strong></td><td>{fields}</td>
            </tr><tr>
                <td><strong>JplaceParser version</strong></td><td>{parser}</td>
            </tr><tr>
                <td><strong>Author</strong></td><td>{author}</td>
            </tr>
        </table>
        """.format(n_plac=len(self.placements),
                   fields=", ".join(self.fields),
                   parser=__version__,
                   author=__author__)
    
    @classmethod
    def fromJplaceFile(cls, jplace: Path) -> JplaceParser:
        """
        Initialize class from path to jplace file
        """
        jplace_obj = cls.getJSONobject(jplace)
        return cls(jplace_obj)

    @staticmethod
    def getJSONobject(json_file: Path) -> dict:
            with open(json_file, 'r') as JSON:
                return json.load(JSON)

    @property
    def meta(self):
        """
        Print metadata
        """
        return self._jplace_obj['metadata']

    @property
    def fields(self): 
        """
        Print data fields
        """
        return self._jplace_obj['fields']
    
    @property
    def placements(self):
        """
        Return placement objects
        """
        return self._jplace_obj['placements']

    def getJplace(self) -> dict:
        return self._jplace_obj

    def writeToFile(self, output_file: Path) -> None:
        """
        Write Jplace object to file
        """
        with open(output_file, 'w') as ofile:
            json.dump(self._jplace_obj, ofile, indent=2)

    def getTreeStr(self, newick=False) -> str:
        """
        Return tree string in original or newick format
        Original format contains branch labels
        in curly brackets. Newick format removes
        these labels.
        """
        if newick:
            return self.newickfyTree(self._jplace_obj['tree'])
        else:
            return self._jplace_obj['tree']
    
    @staticmethod
    def newickfyTree(tree_str: str) -> str:
        """
        Remove branch IDs from jplace tree string
        """
        subs_tree = re.sub("\{(\d+)\}", '', tree_str)
        return subs_tree
        # return next(Phylo.parse(StringIO(subs_tree), 'newick'))

    def getReferenceSequences(self) -> list:
        """
        Get list of reference sequences in the placement tree
        """
        return [c.name for c in self._tree_obj.get_terminals()]
    
    def buildBranchDict(self) -> dict:
        """
        Build dictionary with edge/branch numbers as keys and 
        reference tree leaves as values
        """
        def get_id(s):
            return int(re.search("\{(\d+)\}", s).group(1))
        
        original_tree = self._jplace_obj['tree']
        leaves = self._tree_obj.get_terminals()

        branches = {
            get_id(original_tree[original_tree.find(leaf.name):]): leaf.name
            for leaf in leaves
        }
        return branches

    def _extractPlacementFields(self, pfielddata: list) -> dict:
        fields = self._jplace_obj['fields']
        return {field: pfielddata[i] for i, field in enumerate(fields)}

    def selectBestPlacement(self, placement_object: dict) -> dict:
        """
        Select placement with lowest likelihood
        """
        pdata = [
            self._extractPlacementFields(pfielddata)
            for pfielddata in placement_object['p']
            ]
        lowest_like_placement = sorted(pdata, key=lambda x: x['likelihood'])[0]
        return {'p': lowest_like_placement, 'n': placement_object['n']}

    def selectBestPlacements(self):
        """
        Select placement with lowest likelihood for 
        all placement objects in placements
        """
        best_placements = [
            self.selectBestPlacement(placement)
            for placement in self._jplace_obj['placements']
        ]
        return best_placements

    def computeTreeDiameter(self) -> float:
        """
        Find maximum (pairwise) distance between two tips
        (leaves) in the tree
        """
        root = self._tree_obj.root
        max_distance = 0.0
        tips = self._tree_obj.get_terminals()
        for tip in tips:
            self._tree_obj.root_with_outgroup(tip)
            new_max = max(self._tree_obj.depths().values())
            if new_max > max_distance:
                max_distance = new_max
        self._tree_obj.root_with_outgroup(root)
        return max_distance

    def filterPlacementsByMinimumLWR(self, minimum_lwr: float, outfile: str = None) -> None:
        """
        Filter placements by minimum LWR (from 0 to 1)
        """
        if outfile is None:
            base, ext = os.path.splitext(self._path_to_jplace)
            outfile = f'{base}_min_lwr_{minimum_lwr}{ext}'
        print(f'Filtering placements by minimum LWR: {minimum_lwr}')
        jplace = self.getJSONobject()

        filtered_placement_objs= []
        for placement_object in jplace['placements']:
            filtered_placements = []
            for placement in placement_object['p']:
                edge_num, likelihood, lwr, distal_length, pendant_length = placement
                if lwr >= minimum_lwr:
                    filtered_placements.append(placement)
            if filtered_placements:
                placement_object['p'] = filtered_placements
                filtered_placement_objs.append(placement_object)
        jplace['placements'] = filtered_placement_objs
        with open(outfile, 'w') as ofile:
            json.dump(jplace, ofile, indent=2)

    def filterByMaxPendantToTreeDiameterRatio(self, max_pendant_ratio: float) -> JplaceParser:
        """
        Filter placements by maximum pendant length
        """
        tree_diameter = self.computeTreeDiameter()
        print(f'Filtering placements for tree diameter: {tree_diameter}')
        filtered_jplace = self.getJplace()
        filtered_placement_objs= []
        for placement_object in filtered_jplace['placements']:
            filtered_placements = []
            for placement in placement_object['p']:
                edge_num, likelihood, lwr, distal_length, pendant_length = placement
                if pendant_length / tree_diameter <= max_pendant_ratio:
                    filtered_placements.append(placement)
            if filtered_placements:
                placement_object['p'] = filtered_placements
                filtered_placement_objs.append(placement_object)
        filtered_jplace['placements'] = filtered_placement_objs
        return JplaceParser(filtered_jplace)

    def filterByMaxPendantLength(self, max_pendant_length: float) -> JplaceParser:
        """
        Filter placements by maximum pendant length
        """
        filtered_jplace = self.getJplace()
        filtered_placement_objs= []
        for placement_object in filtered_jplace['placements']:
            filtered_placements = []
            for placement in placement_object['p']:
                edge_num, likelihood, lwr, distal_length, pendant_length = placement
                if pendant_length <= max_pendant_length:
                    filtered_placements.append(placement)
            if filtered_placements:
                placement_object['p'] = filtered_placements
                filtered_placement_objs.append(placement_object)
        filtered_jplace['placements'] = filtered_placement_objs
        return JplaceParser(filtered_jplace)

    def filterByMaxPendantToDistalLengthRatio(self, max_pendant_ratio: float) -> JplaceParser:
        """
        Filter placements by maximum pendant length
        """
        filtered_jplace = self.getJplace()
        filtered_placement_objs= []
        for placement_object in filtered_jplace['placements']:
            filtered_placements = []
            for placement in placement_object['p']:
                edge_num, likelihood, lwr, distal_length, pendant_length = placement
                if pendant_length / distal_length <= max_pendant_ratio:
                    filtered_placements.append(placement)
            if filtered_placements:
                placement_object['p'] = filtered_placements
                filtered_placement_objs.append(placement_object)
        filtered_jplace['placements'] = filtered_placement_objs
        return JplaceParser(filtered_jplace)