#write an option '-clear' to remove intermediate files
$IN_GB_FILE = Convert-Path '.\data\norovirus.gb'
$MIN_PARCER = 0 #nt
$MAX_PARCER = 8000 #nt
$MIN_REMOVE_SIMILAR = '0.5' #percentage e.g. '0.5'=0.5%
$MAX_REMOVE_SIMILAR = '100' #percentage e.g. '100'=100%
$GAPS_IN_ROW_TO_REMOVE = '100' #if a fasta-sequence contains n gaps in a row, it will be removed
$PATH_SCRIPTS = Convert-Path '.\scripts\'
$COORD_FILE = Convert-Path '.\data\norovirus_orf.txt'

$path_dir = (($IN_GB_FILE.split('\'))[0..(($IN_GB_FILE.split('\')).GetUpperBound(0)-1)] -join '\') + '\'

Set-Location ($PATH_SCRIPTS + '..')

Out-Host -InputObject '==========> PARSING'

python ($PATH_SCRIPTS + 'parser_gb.py') -input $IN_GB_FILE -min $MIN_PARCER -max $MAX_PARCER -f 'country,host,collection_date'

Out-Host -InputObject '==========> RESOLVE AMBIGUOUS NT'

python ($PATH_SCRIPTS + 'resolve_ambiguous.py') -input ($IN_GB_FILE.split('.')[0] + '.fasta')

Out-Host -InputObject '==========> SPLIT GENOME'

python ($PATH_SCRIPTS + 'split_genome.py') -input ($IN_GB_FILE.split('.')[0]  + '_less_amb.fasta') -coord $COORD_FILE

$coord_arr = ((Get-Content -path $COORD_FILE -head 1).split(','))
$coord_arr = $coord_arr[1..($coord_arr.Length-1)]

Out-Host -InputObject '==========> ALIGNMENTS'

foreach ($frame in $coord_arr)
{
	python ($PATH_SCRIPTS + 'trans_alignment.py') -input ($IN_GB_FILE.split('.')[0] + '_less_amb_' + $frame + '.fasta')
	trimal -gt 0.8 -in ($IN_GB_FILE.split('.')[0] + '_less_amb_' + $frame + '_tr_aln_rt.fasta') -out ($IN_GB_FILE.split('.')[0] + '_less_amb_' + $frame + '_tr_aln_rt_gap-0.2.fasta')
}
$join_list = @()
foreach ($frame in $coord_arr)
{
	$join_list += ($IN_GB_FILE.split('.')[0] + '_less_amb_' + $frame + '_tr_aln_rt_gap-0.2.fasta')
}
$join_list = $join_list -join ','

Out-Host -InputObject '==========> JOIN ALIGNMENTS'

python ($PATH_SCRIPTS + 'join_al.py') -input_list $join_list -out_name ($IN_GB_FILE.split('.')[0] + '_full_aln.fasta')

Out-Host -InputObject '==========> REMOVE GAPS'

python ($PATH_SCRIPTS + 'gap_in_row.py') -input ($IN_GB_FILE.split('.')[0] + '_full_aln.fasta') -gap_count $GAPS_IN_ROW_TO_REMOVE

Out-Host -InputObject '==========> REMOVE SIMILAR'

python ($PATH_SCRIPTS + 'remove_similar.py') -input ($IN_GB_FILE.split('.')[0] + '_full_aln_' + $GAPS_IN_ROW_TO_REMOVE + 'gp.fasta') -min $MIN_REMOVE_SIMILAR -max $MAX_REMOVE_SIMILAR
#mkdir $path_dir + 'genotyped'
#cp ($IN_GB_FILE.split('.')[0] + '_full_aln_' + $GAPS_IN_ROW_TO_REMOVE + 'gp_' + $MIN_REMOVE_SIMILAR + '.fasta') ($path_dir + 'genotyped')
#python ($PATH_SCRIPTS + 'genotyping.py') 
if ($args[0] -eq '-clear')
{
	rm ($path_dir + 'norovirus_less_amb*')
	rm ($path_dir + 'norovirus_full_aln.fasta')
	rm ($path_dir + 'norovirus.fasta')
	rm ($path_dir + 'norovirus_slices.fasta')
	rm ($path_dir + 'local_db*')
	rm ($path_dir + 'blast.out')
}