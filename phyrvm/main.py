import logging
import os
import sys
import subprocess
from multiprocessing import Process
from tqdm import tqdm
from shutil import copy
from phyrvm.config.config import LOG_TASK
from phyrvm.rvm.trim_contamination import Trim_contamination
from phyrvm.rvm.summary import Summary
from phyrvm.rvm.rmdup import Rmdup
from phyrvm.rvm.calculate_RPM import Calculate_RPM
from .phylogenetic_analysis import Phylogenetic_analysis
from phyrvm.biolib.common import (check_dependencies,
                                  check_file_exists,
                                  make_sure_path_exists)
from phyrvm.biolib.exceptions import PhyRvmException
from .contigs_filter import Contigs_filter

class OptionsParser(object):

    def __init__(self, version):
        self.logger = logging.getLogger('timestamp')
        self.warnings = logging.getLogger('warnings')
        self.version = version
        self.genomes_to_process = None


    def contigs_filter(self,options):

        self.logger.log(LOG_TASK,'START Primary Screen PART...')

        contigs_filter_path = os.path.join(options.out_dir,'contigs_filter')
        res_path = os.path.join(options.out_dir,"results")
        contigs_filter_fasta = f"{res_path}/Primary_screen_res.fasta"
        contigs_filter_fasta_keep_dup = f"{res_path}/Primary_screen_keep_dup_res.fasta"
        check_file_exists(options.i)
        make_sure_path_exists(contigs_filter_path)
        make_sure_path_exists(os.path.join(options.out_dir,'results'))

        if options.i2 is not None:
            check_file_exists(options.i2)
            reads_list = [options.i,options.i2]
        else:
            reads_list = [options.i]

        contigs_filter_item = Contigs_filter(reads_list,contigs_filter_path,options.threads)
        rm_rRNA_file,filename = contigs_filter_item.qc_with_rm_rRNA()

        accession_tax_VirusesFlitered_file = contigs_filter_item.run(options,rm_rRNA_file,filename)

        if len(reads_list) == 2:
            os.remove(os.path.join(contigs_filter_path,"step1_QC_1.fq.gz"))
            os.remove(os.path.join(contigs_filter_path,"step1_QC_2.fq.gz"))
            os.remove(os.path.join(contigs_filter_path,"step2_QC_1.fq"))
            os.remove(os.path.join(contigs_filter_path,"step2_QC_2.fq"))
            os.remove(os.path.join(contigs_filter_path,"step3_QC_cdhit_1.fq"))
            os.remove(os.path.join(contigs_filter_path,"step3_QC_cdhit_1.fq.clstr"))
            os.remove(os.path.join(contigs_filter_path,"step3_QC_cdhit_1.fq2.clstr"))
            os.remove(os.path.join(contigs_filter_path,"step3_QC_cdhit_2.fq"))
        else:
            os.remove(os.path.join(contigs_filter_path,"step1_QC.fq.gz"))
            os.remove(os.path.join(contigs_filter_path,"step2_QC.fq"))
            os.remove(os.path.join(contigs_filter_path,"step3_QC_cdhit.fq"))

        if not os.path.isfile(accession_tax_VirusesFlitered_file):
            self.logger.log(LOG_TASK,'END Primary Screen PART...')
            sys.exit()
        elif os.path.getsize(accession_tax_VirusesFlitered_file) < 1:
            self.logger.log(LOG_TASK,'END Primary Screen PART...')
            sys.exit()

        if options.no_trim_contamination is False:
            self.logger.info('[contigs_filter] Cut the sequence contamination at both ends of contigs')
            output_nt_fasta = Trim_contamination(contigs_filter_path,options.threads).run(accession_tax_VirusesFlitered_file)
            res_file = output_nt_fasta
        else:
            res_file = accession_tax_VirusesFlitered_file

        if not os.path.isfile(res_file):
            self.logger.log(LOG_TASK,'END Primary Screen PART...')
            sys.exit()
        elif os.path.getsize(res_file) < 1:
            self.logger.log(LOG_TASK,'END Primary Screen PART...')
            sys.exit()

        rmdup_res_fas = Rmdup(options.subparser_name,contigs_filter_path,options.threads).run(res_file)
        copy(os.path.join(contigs_filter_path,'step11_rmdup_res.fas'),contigs_filter_fasta)
        rmdup_res_fas_keep_dup = res_file
        copy(rmdup_res_fas_keep_dup,contigs_filter_fasta_keep_dup)
        
        Summary(options,contigs_filter_path,"Primary_screen_res.tsv").run()
        Summary(options,contigs_filter_path,"Primary_screen_keep_dup_res.tsv").run()

        Calculate_RPM(contigs_filter_path,
                    options.threads,"Primary_screen_res.tsv").run(rmdup_res_fas,rm_rRNA_file)
        Calculate_RPM(contigs_filter_path,
            options.threads,"Primary_screen_keep_dup_res.tsv").run(rmdup_res_fas_keep_dup,rm_rRNA_file)
        
        self.logger.log(LOG_TASK,'END Primary Screen PART...')

        return rmdup_res_fas



    def phylogenetic_analysis(self,options):
        self.logger.log(LOG_TASK,'START Contigs Classify  PART...')
        if os.path.getsize(options.classify_i) < 1:
            self.logger.warn(f'There is no sequence in the file.')
            return 0
        contigs_classify_path = os.path.join(options.out_dir,'phylogenetic_analysis')
        make_sure_path_exists(contigs_classify_path)
        make_sure_path_exists(os.path.join(options.out_dir,'results'))
        check_file_exists(options.classify_i)

        if options.keep_dup is False:
            rmdup_res_fas = Rmdup(options.subparser_name,
                        contigs_classify_path,options.threads).run(options.classify_i)
        else:
            rmdup_res_fas = options.classify_i

        classify_item = Phylogenetic_analysis(rmdup_res_fas,
                                         contigs_classify_path,options)
        classify_item.run()
        self.logger.log(LOG_TASK,'END Contigs Classify  PART...')


    def parse_options(self, options):
        if sys.version_info.major < 3:
            raise PhyRvmException('Python 2 is no longer supported.')

        if hasattr(options, 'threads') and options.threads < 1:
            self.logger.warning(
                'You cannot use less than 1 CPU, defaulting to 1.')
            options.cpus = 1

        if options.subparser_name == 'end_to_end':
            check_dependencies(['bbduk.sh',
                    'diamond','seqkit','bowtie2','cd-hit','taxonkit',
                    'megahit','makeblastdb','blastn',
                    'blastp','mafft','trimal','phyml','pplacer','guppy'])

            contigs_file = self.contigs_filter(options)
            options.classify_i = contigs_file
            self.phylogenetic_analysis(options)

        elif options.subparser_name == 'contigs_filter':
            check_dependencies(['bbduk.sh','cd-hit',
                    'diamond','seqkit','bowtie2','taxonkit',
                    'megahit','makeblastdb','blastn',
                    'blastp'])

            self.contigs_filter(options)
        elif options.subparser_name == 'phylogenetic_analysis':
            check_dependencies(['bbduk.sh','cd-hit',
                    'diamond','seqkit','makeblastdb','blastn','taxonkit',
                    'blastp','mafft','trimal','phyml','pplacer','guppy'])

            self.phylogenetic_analysis(options)

        elif options.subparser_name == 'install':
            parent_path = os.path.dirname(os.path.realpath(__file__))
            args = f'pip install -r {parent_path}/requirements.txt'
            subprocess.call(args)
        else:
            self.logger.error('Unknown phyrvm command: "' +
                              options.subparser_name + '"\n')
            sys.exit(1)

        return 0

