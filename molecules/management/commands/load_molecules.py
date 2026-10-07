import csv
from django.core.management.base import BaseCommand
from molecules.models import Molecule

class Command(BaseCommand):
    help = 'Imports cleaned molecular data from a CSV file into the database'

    def add_arguments(self, parser):
        parser.add_argument('csv_file', type=str, help='Path to the cleaned_molecules.csv file')

    def handle(self, *args, **options):
        csv_file_path = options['csv_file']
        
        self.stdout.write(self.style.SUCCESS(f"Reading data from {csv_file_path}..."))
        
        molecules_to_create = []
        
        try:
            with open(csv_file_path, mode='r', encoding='utf-8') as file:
                reader = csv.DictReader(file)
                
                for row in reader:
                    # Map CSV columns straight to your Django model fields
                    mol = Molecule(
                        molecule_id=row['molecule_id'],
                        smiles=row['smiles'],
                        molecular_weight=float(row['molecular_weight']),
                        log_p=float(row['log_p']),
                        tpsa=float(row['tpsa']),
                        h_bond_donors=int(row['h_bond_donors']),
                        h_bond_acceptors=int(row['h_bond_acceptors']),
                        rotatable_bonds=int(row['rotatable_bonds']),
                        ring_count=int(row['ring_count']),
                        fingerprint=row['fingerprint']
                    )
                    molecules_to_create.append(mol)
                    
            # Break up the insert into smaller chunks so the Free Tier database doesn't timeout
            self.stdout.write(f"Inserting {len(molecules_to_create)} records into database in batches...")
            Molecule.objects.bulk_create(molecules_to_create, batch_size=500, ignore_conflicts=True)
            self.stdout.write(self.style.SUCCESS(f"Successfully loaded molecules."))

            
        except FileNotFoundError:
            self.stdout.write(self.style.ERROR(f"File not found at: {csv_file_path}"))
