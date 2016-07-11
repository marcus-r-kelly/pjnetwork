from django.db import connection
from django.core.management.base import BaseCommand

from lib.markutils import absoverlap
from network.models import Entrez


symbols = ['C9orf142', 'MRPL50', 'ZNF784', 'USP44', 'TGS1', 'ND5', 'PRSS35', 'ACADS', 'ACADSB', 'ACADVL', 'ACAT1', 'ACAT2', 'ASIC2', 'ASIC1', 'ACHE', 'ACLY', 'ACO1', 'ACR', 'ACO2', 'ACOX1', 'ACP1', 'ACP2', 'ACP5', 'ACPP', 'ACRV1', 'ACTA1', 'ACTA2', 'ACTB', 'ACTC1', 'ACTG1', 'ACTG2', 'ACTN4', 'ACTL6A', 'ACTN1', 'ACTN2', 'ACTN3', 'ACVR1', 'ACVR1B', 'ACVR2A', 'ACVR2B', 'ACVRL1', 'ACY1', 'ACYP1', 'ACYP2', 'ADA', 'ADAM8', 'ADAM10', 'ADAR', 'ADARB1', 'ADARB2', 'ADCY1', 'ADCY2', 'ADCY3', 'ADCY5', 'ADCY6', 'ADCY7', 'ADCY8', 'ADCY9', 'ADCYAP1', 'ADCYAP1R1', 'ADD1', 'ADD2', 'ADD3', 'PLIN2', 'ADH1A', 'ADH1B', 'ADH1C', 'ADH4', 'ADH5', 'ADH6', 'ADH7', 'ADK', 'ADM', 'ADORA1', 'ADORA2A', 'ADORA2B', 'ADORA3', 'ADPRH', 'PARP1', 'PARP4', 'ADRA1D', 'ADRA1B', 'ADRA1A', 'ADRA2A', 'ADRA2B', 'ADRA2C', 'ADRB1', 'ADRB2', 'ADRB3', 'ADRBK1', 'ADRBK2', 'ADSL', 'ADSS', 'AP2A1', 'AP2A2', 'AP1B1', 'AP2B1', 'AP1G1', 'AEBP1', 'AES', 'CRISP1', 'AFM', 'AFP', 'AGA', 'ACAN', 'AGER', 'AGL', 'AGRP', 'JAG1', 'AGT', 'AGTR1', 'AGTR2', 'APLNR', 'AGXT', 'NR0B1', 'AHCY', 'AHR', 'AHSG', 'AIF1', 'AIM1', 'AK1', 'AK2', 'AK4', 'AKT1', 'AKT2', 'ALAD', 'ALAS1', 'ALAS2', 'ALB', 'ALCAM', 'kenyer']


class Command(BaseCommand):

    import pyexcel as pe
    help = 'our help string comes here'
    
    def handle(self, *args, **options):

        absoverlap(symbols)
