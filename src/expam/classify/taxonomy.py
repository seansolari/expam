import os
import re
import time
import requests

from expam.classify import _STEP, ENTREZ_API, ENTREZ_EMAIL, ENTREZ_TOOL
from expam.database import FileLocationConfig
from expam.database.config import validate_taxonomy_files
from expam.utils import yield_csv


class TaxonomyNCBI:
    def __init__(self, file_config: FileLocationConfig) -> None:
        self.config: FileLocationConfig = file_config

    def find_downloaded_taxonomy(self):
        if not validate_taxonomy_files(self.config):
            raise OSError("Taxonomy files not located!")

    def load_taxonomy_map(self, convert_to_name=True):
        self.find_downloaded_taxonomy()

        # Create map from scientific name --> (taxid, rank).
        taxon_data = {}
        for data in yield_csv(self.config.taxon_rank):
            taxon_data[",".join(data[0:-2])] = tuple(data[-2:])

        # Create map from tax_id --> lineage (tuple).
        tax_id_to_lineage = {}
        for data in yield_csv(self.config.taxid_lineage):
            tax_id_to_lineage[data[0]] = tuple(data[1:])

        if not convert_to_name:
            return tax_id_to_lineage, taxon_data

        # Create map from name --> lineage (tuple).
        name_to_lineage = {}
        for data in yield_csv(self.config.accession_id):
            name_to_lineage[data[0]] = tax_id_to_lineage[data[2]]

        return name_to_lineage, taxon_data

    def load_sequence_map(self):
        return list(yield_csv(self.config.accession_id))

    def load_taxid_lineage_map(self):
        if os.path.exists(self.config.taxid_lineage):
            return list(yield_csv(self.config.taxid_lineage))
        else:
            return []

    def load_rank_map(self):
        name_to_rank = {}
        
        if os.path.exists(self.config.taxon_rank):
            for data in yield_csv(self.config.taxon_rank):
                if len(data) > 1:
                    name_to_rank[data[0]] = ",".join(data[1:])
        
        return name_to_rank

    def accession_to_taxonomy(self):
        """
        Map accession IDs to taxonomic labels.

        :sequence_ids: List - (sequence_id, accession_id, taxon_id)
        :taxon_ranks: List - (taxon_id, taxonomic ranks)
        :taxa_to_rank: Dict - taxon --> rank
        """

        def tuples_to_disk(lst):
            return "\n".join([",".join(item) for item in lst])

        def dict_to_disk(dct):
            return "\n".join([",".join((key, value)) for key, value in dct.items()])

        sequence_ids = self.load_sequence_map()
        taxon_ranks = self.load_taxid_lineage_map()
        taxa_to_rank = self.load_rank_map()

        # Collect taxon ids for unknown organisms.
        accessions_to_be_mapped = []
        taxa_to_be_collected = []

        for (sequence_id, accession_id, taxon_id) in sequence_ids:
            if taxon_id == "None":
                accessions_to_be_mapped.append(accession_id)
            else:
                taxa_to_be_collected.append(taxon_id)

        if accessions_to_be_mapped:
            requestor = EntrezRequestor()
            accession_to_tax = requestor.request_tax_ids("nuccore", accessions_to_be_mapped)

            print("Received %d response(s) for ESummary TaxID request!"
                % len(accession_to_tax))

            for i in range(len(sequence_ids)):
                if sequence_ids[i][1] in accession_to_tax:
                    sequence_ids[i][2] = accession_to_tax[sequence_ids[i][1]]
                    taxa_to_be_collected.append(sequence_ids[i][2])

        # Collect taxonomic lineages for taxa.
        current_taxa = {taxa[0] for taxa in taxon_ranks}
        taxa_to_be_collected = {  # Set so that we collect unique values.
            taxon_id
            for taxon_id in taxa_to_be_collected
            if taxon_id not in current_taxa
        }

        if taxa_to_be_collected:
            taxid_to_taxon, taxon_to_rank = requestor.request_labels("taxonomy", "xml", list(taxa_to_be_collected))

            print("Received %d response(s) for EFetch Taxon request!"
                % len(taxid_to_taxon))

            taxon_ranks.extend(taxid_to_taxon)
            taxa_to_rank.update(taxon_to_rank)

        # Save update maps to disk.
        with open(self.config.accession_id, "w") as f:
            f.write(tuples_to_disk(sequence_ids))

        with open(self.config.taxid_lineage, "w") as f:
            f.write(tuples_to_disk(taxon_ranks))

        print("Taxonomic lineages written to %s!" % self.config.taxid_lineage)

        with open(self.config.taxon_rank, "w") as f:
            f.write(dict_to_disk(taxa_to_rank))

        print("Taxonomic ranks written to %s!" % self.config.taxon_rank)


class EntrezRequestor:
    def __init__(self, entrez_tool: str = None, entrez_email: str = None, api_key: str = None) -> None:
        self.entrez_tool = ENTREZ_TOOL if entrez_tool is None else entrez_tool
        self.entrez_email = ENTREZ_EMAIL if entrez_email is None else entrez_email
        self.api_key = ENTREZ_API if api_key is None else api_key

    def request_tax_ids(self, db, id_list):
        POST_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        taxids, upto = {}, 0

        PARAMS = {
            "tool": self.entrez_tool,
            "email": self.entrez_email,
            "api_key": self.api_key,
            "db": db,
        }

        while upto < len(id_list):
            next_requests = id_list[upto:upto + _STEP]

            print("Posting %d UIDs to NCBI Entrez %s."
                % (len(next_requests), db))

            PARAMS["id"] = ",".join(next_requests)
            esummary_request = requests.post(
                url=POST_URL,
                data=PARAMS
            )

            # Parse TaxonIDs from raw results.
            accn_id = tax_id = None

            for line in esummary_request.text.split("\n"):
                if "<DocSum>" in line:  # Start new record.
                    accn_id = tax_id = None

                elif 'AccessionVersion' in line:
                    accn_id = self.parse_id(line)

                elif 'TaxId' in line:
                    tax_id = self.parse_tax_id(line)

                elif "</DocSum>" in line:  # Only save complete records.
                    if accn_id is not None and tax_id is not None:
                        taxids[accn_id] = tax_id

            upto += _STEP

            time.sleep(1.0)  # Allow server time to breath.

        return taxids

    @staticmethod
    def parse_id(string):
        new_id = re.findall(r'\<Item Name\="AccessionVersion" Type\="String"\>(.*?)\<\/Item\>', string)

        if not new_id:
            raise ValueError("No taxids found!")
        else:
            return new_id[0]

    @staticmethod
    def parse_tax_id(string):
        taxids = re.findall(r'\<Item Name\="TaxId" Type\="Integer"\>(.*?)\<\/Item\>', string)

        if not taxids:
            raise ValueError("No taxids found!")
        else:
            return taxids[0]

    def request_labels(self, db, retmode, id_list):
        POST_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        taxon_labels, ranks, upto = [], {}, 0

        PARAMS = {
            "tool": self.entrez_tool,
            "email": self.entrez_email,
            "api_key": self.api_key,
            "db": db,
            "retmode": retmode,
        }

        while upto < len(id_list):
            next_requests = id_list[upto:upto + _STEP]

            print("Posting %d UIDs to NCBI Entrez %s."
                % (len(next_requests), db))

            PARAMS["id"] = ",".join(next_requests)
            efetch_request = requests.post(
                url=POST_URL,
                data=PARAMS
            )

            # Parse taxonomic labels from raw results.
            tax_id = taxa = rank = name = None
            sub_tax_id = sub_rank = sub_name = None
            collect = True

            for line in efetch_request.text.split("\n"):
                if "<Taxon>" in line and collect:  # Start new record.
                    tax_id = taxa = name = rank = None

                elif "<Taxon>" in line and not collect:
                    sub_tax_id = sub_rank = sub_name = None

                elif "<TaxId>" in line and collect:
                    tax_id = self.parse_single_tax_id(line)

                elif "<TaxId>" in line and not collect:
                    sub_tax_id = self.parse_single_tax_id(line)

                elif "<Lineage>" in line and collect:
                    taxa = self.parse_lineage(line).replace(",", "")

                elif "<Rank>" in line and not collect:
                    if sub_name == "cellular organisms":
                        sub_rank = "top"

                    else:
                        sub_rank = self.parse_rank(line)

                elif "<Rank>" in line and collect:
                    if name == "cellular organisms":
                        rank = "top"

                    else:
                        rank = self.parse_rank(line)

                elif "<ScientificName>" in line and not collect:
                    sub_name = self.parse_name(line)

                elif "<ScientificName>" in line and collect:
                    name = self.parse_name(line)

                elif "<LineageEx>" in line:
                    collect = False

                elif "</Taxon>" in line and not collect:
                    ranks[sub_name] = sub_tax_id + "," + sub_rank

                elif "</LineageEx>" in line:
                    collect = True

                elif "</Taxon>" in line and collect:
                    if tax_id is not None and taxa is not None:
                        lineage = taxa.strip().split("; ")
                        if name not in lineage and name is not None:
                            lineage += [name]

                        taxon_labels.append([tax_id, ",".join(lineage)])

                    ranks[name] = tax_id + ',' + rank

            upto += _STEP

            time.sleep(1.0)  # Allow server time to breath.

        return taxon_labels, ranks

    @staticmethod
    def parse_single_tax_id(string):
        taxids = re.findall(r'\<TaxId\>(.*?)\<\/TaxId\>', string)

        if not taxids:
            raise ValueError("No taxids found!")
        else:
            return taxids[0]

    @staticmethod
    def parse_lineage(string):
        lineage = re.findall(r'\<Lineage\>(.*?)\<\/Lineage\>', string)

        if not lineage:
            raise ValueError("No lineages found!")
        else:
            return lineage[0]

    @staticmethod
    def parse_rank(string):
        rank = re.findall(r'\<Rank\>(.*?)\<\/Rank\>', string)

        if not rank:
            raise ValueError("No rank found!")
        else:
            return rank[0]

    @staticmethod
    def parse_name(string):
        name = re.findall(r'\<ScientificName\>(.*?)\<\/ScientificName\>', string)

        if not name:
            raise ValueError("No rank found!")
        else:
            name = re.sub(r"[,]", "", name[0])
            return name



