from django.db import models

# Create your models here.

class Molecule(models.Model):
    molecule_id = models.CharField(max_length=100, unique=True, help_text="Molecule ID from ANPDB")
    smiles = models.TextField(help_text="SMILES string representation")
    
    # Physicochemical properties
    molecular_weight = models.FloatField()
    log_p = models.FloatField(help_text="Lipophilicity")
    tpsa = models.FloatField(help_text="Polarity / Total Polar Surface Area")
    h_bond_donors = models.IntegerField()
    h_bond_acceptors = models.IntegerField()
    rotatable_bonds = models.IntegerField()
    ring_count = models.IntegerField()
    
    # Store the 2048-bit fingerprint as text for direct string comparison
    fingerprint = models.TextField(help_text="Bit string representation of Morgan FP")
    
    created_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return f"{self.molecule_id} ({self.smiles[:15]}...)"
