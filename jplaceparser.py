import os
import re 
import json 
from io import StringIO

from Bio import Phylo 


class JplaceParser():
    """
    Methods to parse jplace files, as specified in 
    https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0031009
    """
    def __init__(self, path_to_jplace: str) -> None:
        self._path_to_jplace = path_to_jplace
        self.jplace = self.getJSONobject()
        self._tree_obj = next(Phylo.parse(StringIO(self.newickfyTree(self.jplace['tree'])), 'newick'))

    def getJSONobject(self) -> dict:
        with open(self._path_to_jplace, 'r') as JSON:
            return json.load(JSON)
    
    @property
    def meta(self):
        """
        Print metadata
        """
        return self.jplace['metadata']

    @property
    def fields(self): 
        """
        Print data fields
        """
        return self.jplace['fields']
    
    @property
    def placements(self):
        """
        Return placement objects
        """
        return self.jplace['placements']

    def getTreeStr(self, newick=False) -> str:
        """
        Return tree string in original or newick format
        Original format contains branch labels
        in curly brackets. Newick format removes
        these labels.
        """
        if newick:
            return self.newickfyTree(self.jplace['tree'])
        else:
            return self.jplace['tree']
    
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
        
        original_tree = self.jplace['tree']
        leaves = self._tree_obj.get_terminals()

        branches = {
            get_id(original_tree[original_tree.find(leaf.name):]): leaf.name
            for leaf in leaves
        }
        return branches

    def extractPlacementFields(self, pfielddata: list) -> dict:
        """
        Get dict with placement field values from list of values
        """
        fields = self.jplace['fields']
        return {field: pfielddata[i] for i, field in enumerate(fields)}

    def selectBestPlacement(self, placement_object: dict) -> dict:
        """
        Select placement with lowest likelihood
        """
        pdata = [
            self.extractPlacementFields(pfielddata)
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
            for placement in self.jplace['placements']
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

    def filterPlacementsByMaxPendantToTreeDiameterRatio(self, max_pendant_ratio: float,
                                                        outfile: str = None) -> None:
        """
        Filter placements by maximum pendant length
        """
        if outfile is None:
            base, ext = os.path.splitext(self._path_to_jplace)
            outfile = f'{base}_max_pendant_diameter_ratio_{max_pendant_ratio}{ext}'
        tree_diameter = self.computeTreeDiameter()
        print(f'Filtering placements for tree diameter: {tree_diameter}')
        jplace = self.getJSONobject()

        filtered_placement_objs= []
        for placement_object in jplace['placements']:
            filtered_placements = []
            for placement in placement_object['p']:
                edge_num, likelihood, lwr, distal_length, pendant_length = placement
                if pendant_length / tree_diameter <= max_pendant_ratio:
                    filtered_placements.append(placement)
            if filtered_placements:
                placement_object['p'] = filtered_placements
                filtered_placement_objs.append(placement_object)
        jplace['placements'] = filtered_placement_objs
        
        with open(outfile, 'w') as ofile:
            json.dump(jplace, ofile, indent=2)

    def filterPlacementsByMaxPendantLength(self, max_pendant_length: float, outfile: str = None) -> None:
        """
        Filter placements by maximum pendant length
        """
        if outfile is None:
            base, ext = os.path.splitext(self._path_to_jplace)
            outfile = f'{base}_max_pendant_{max_pendant_length}{ext}'
        jplace = self.getJSONobject()

        filtered_placement_objs= []
        for placement_object in jplace['placements']:
            filtered_placements = []
            for placement in placement_object['p']:
                edge_num, likelihood, lwr, distal_length, pendant_length = placement
                if pendant_length <= max_pendant_length:
                    filtered_placements.append(placement)
            if filtered_placements:
                placement_object['p'] = filtered_placements
                filtered_placement_objs.append(placement_object)
        jplace['placements'] = filtered_placement_objs
        
        with open(outfile, 'w') as ofile:
            json.dump(jplace, ofile, indent=2)

    def filterPlacementsByMaxPendantToDistalLengthRatio(self, max_pendant_ratio: float,
                                                        outfile: str = None) -> None:
        """
        Filter placements by maximum pendant length
        """
        if outfile is None:
            base, ext = os.path.splitext(self._path_to_jplace)
            outfile = f'{base}_max_pendant_distal_ratio_{max_pendant_ratio}{ext}'
        jplace = self.getJSONobject()

        filtered_placement_objs= []
        for placement_object in jplace['placements']:
            filtered_placements = []
            for placement in placement_object['p']:
                edge_num, likelihood, lwr, distal_length, pendant_length = placement
                if pendant_length / distal_length <= max_pendant_ratio:
                    filtered_placements.append(placement)
            if filtered_placements:
                placement_object['p'] = filtered_placements
                filtered_placement_objs.append(placement_object)
        jplace['placements'] = filtered_placement_objs
        
        with open(outfile, 'w') as ofile:
            json.dump(jplace, ofile, indent=2)