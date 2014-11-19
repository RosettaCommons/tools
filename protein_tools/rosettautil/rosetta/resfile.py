class Residue_File():
    """A class that will handle residue files"""

    def __init__(self, resfile):
        self.res_file_dict = {}
        self._parse(resfile)

    def _parse(self, resfile):
        breaker = False
        self.entity = 1
        resfile_lines = open(resfile).readlines()
        if not resfile_lines:
            raise Exception("RESFILE is empty")
        for res_line in resfile_lines:
            if res_line.split()[0].lower() == "start":
                breaker = True
                continue
            if breaker:
                res_number = res_line.split()[0]
                chain = res_line.split()[1]
                desing_param = res_line.split()[2]
                self.res_file_dict[self.entity] = {"chain": chain, "residue": int(res_number), "design": desing_param}
                self.entity += 1
        if self.entity == 1:
            raise Exception("Your resfile did not parse correctly, is there anything in it after start?")

    def get_designed_entities(self):
        """A method that returns a list of tuples containing the chain and residue pairings of designed entities\
        that where considered for design specified by the resfile"""
        designed_entities = []
        for i in self.res_file_dict:
            designed_entities.append((self.res_file_dict[i]["chain"], self.res_file_dict[i]["residue"]))
        return designed_entities

    def number_designed_entities(self):
        """Returns the number of entities in the resfile"""
        return self.entity
