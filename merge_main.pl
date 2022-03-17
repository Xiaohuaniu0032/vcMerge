use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use FindBin qw/$Bin/;

my ($report_dir,$outdir) = @ARGV;

# SARS_CoV_2_variantCaller_out.1613/
# variantCaller_out.1629/


# first get all variantCaller_out.* dir
my $vc_dirs_aref = &get_vc_dir_by_run_time($report_dir);

# for each variantCaller_out.*, creat a new dir
for my $dir (@{$vc_dirs_aref}){
	my $vc_dir = "$report_dir/plugin_out/$dir";
	#print "$vc_dir\n";
	# get all TSVC_variants.vcf
	#my @vcf = glob "$vc_dir/*/TSVC_variants.vcf";

	my @vcf = glob "$vc_dir/*/alleles.xls"; # alleles_<barcode_xxxx>.xls -> alleles.xls
	# if this vcf is from sars-cov-2 sample?
	# vcf contain "2019-nCoV"
	my @cov_vcf; # SARS vcf
	for my $vcf (@vcf){
		#print "$vcf\n";
		my $if_cov2_vcf = &check_if_cov2_vcf($vcf);
		if ($if_cov2_vcf eq "YES"){
			print "[SARS-CoV-2 VCF]: $vcf\n";
			push @cov_vcf, $vcf;
		}else{
			#print "[NOT SARS-CoV-2 VCF]: $vcf\n";
			next;
		}
	}

	my %sample_vars; # 记录特定样本的变异
	my %all_vars;    # 记录这个run下所有样本的变异
	my %sample;      # 样本名称(barcode编号)
	for my $vcf (@cov_vcf){
		my $base_dir = dirname($vcf);
		my $barcode = basename($base_dir); # IonXpress_001 / IonCode_0109
		$sample{$barcode} = 1; # barcode [IonXpress_001]

		open VCF, "$vcf" or die;
		<VCF>; # skip header line
		while (<VCF>){
			chomp;
			#next if (/^\#/); #skip # line
			my @arr = split /\t/;
			if ($arr[0] ne "2019-nCoV"){
				next; # only keep 2019-nCoV chrom
			}

			next if ($arr[4] ne "Homozygous" and $arr[4] ne "Heterozygous");
			next if ($arr[5] ne "-"); # Filter

			my $freq = $arr[6];
			my $var = "$arr[0]\:$arr[1]\:$arr[2]\:$arr[3]"; # 2019-nCoV:210:G:T
			$sample_vars{$barcode}{$var} = $freq;
			my $pos = $arr[1];
			push @{$all_vars{$pos}}, $var; # 1) one pos may has diff variants; 2) may contain dup vars
		}
		close VCF;
	}
			

	# outdir is /results/analysis/output/Home/2019-nCoV-map2hg19-exon-virus_241/plugin_out/variantCaller_out.1535
	#my $plugin_number = basename($outdir);
	my $outfile = "$outdir/$dir\.TSVC_variants.merged.vcf.xls";
	print "[Merge VCF is]: $outfile\n";
	# Chrom/Position/Ref/Variant/IonXpress_001.Freq/IonXpress_002.Freq/.../
	# 2019-nCoV/210/G/T/1/0.98/.../
	
	open O, ">$outfile" or die;
	print O "Chrom\tPosition\tRef\tVariant";

	my @barcode = keys %sample;
	my @barcode_sort = sort {$a cmp $b} @barcode;

	for my $s (@barcode_sort){
		print O "\t$s";
	}
	print O "\n";

	# sort variant by pos
	foreach my $pos (sort { $a <=> $b } keys %all_vars){
		my $var_aref = $all_vars{$pos}; # one pos may have multi var type
		my %uniq_var;
		for my $var (@{$var_aref}){
			$uniq_var{$var} = 1;
		}

		for my $var (keys %uniq_var){ # 2019-nCoV:210:G:T
			my @var = split /\:/, $var;
			print O "$var[0]\t$var[1]\t$var[2]\t$var[3]";
			for my $s (@barcode_sort){
				my $freq;
				if (exists $sample_vars{$s}{$var}){
					$freq = $sample_vars{$s}{$var};
				}else{
					$freq = "NA";
				}
				print O "\t$freq";
			}
			print O "\n";
		}
	}
	close O;
}



sub check_if_cov2_vcf{
	my $vcf = $_[0];
	my $res;
	open IN, "$vcf" or die;
	<IN>; # skip first line;
	my $first_line = <IN>;
	close IN;

	my $chrom = (split /\t/, $first_line)[0];
	if ($chrom eq "2019-nCoV"){
		$res = "YES";
	}else{
		$res = "NO";
	}

	return($res);
}



sub get_vc_dir_by_run_time{
	my ($dir) = @_;
	my @startplugin_json = glob "$dir/plugin_out/*/startplugin.json"; # /results/analysis/output/Home/2019-nCoV-map2hg19-exon-virus_241/plugin_out/variantCaller_out.1257/startplugin.json
	my @vc_dirs;
	for my $json (@startplugin_json){
		my $basedir = dirname($json);
		my $vc_name = basename($basedir);
		if ($vc_name =~ /variantCaller/){
			next if ($vc_name =~ /variantCallerMerge/);
			# variantCaller_out.1257/
			# SARS_CoV_2_variantCaller_out.1529/
			push @vc_dirs, $vc_name;
		}
	}

	return(\@vc_dirs);
}
