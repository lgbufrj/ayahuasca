#!/usr/bin/perl -w 

########
# Este programa pega um arquivo multifasta e um arquivo contendo o accession number das sequencias desejadas
# no arquivo original. Ele passa todo o arquivo original e gera um novo arquivo contendo apenas
# as sequencias cujo accession foi dado no segundo arquivo
########

if (!$ARGV[1]) { 
  print "Usage: ./get_seq_by_acc.pl <multiFASTA.file> <file containing the accession number of sequences desired (one by line)>\n"; 
  exit; 
}

# Obtem o arquivo de numeros de sequencias a serem obtidas
open (ACCS, $ARGV[1]);
my @acc = ();
while (<ACCS>) {
  chomp;
  push (@acc, $_);
}
close ACCS;

my $printing = 0;

# Abre o arquivo de sequencias
open (MULTIFASTA, $ARGV[0]);
while (<MULTIFASTA>) {
  chomp ($_);
  
  # Para cada nova sequencia, incrementa o contador de sequencia
  if ($_ =~ /^>/) { 
    $printing = 0;
    foreach my $an (@acc) {
      if ($_ =~ $an) { $printing = 1; last; }
    }
  }

  # Enquanto o contador de sequencia for igual ao contador de numero 
  # de sequencia desejado, imprime as linhas.
  if ($printing) {
      print "$_\n";
  } 
}
