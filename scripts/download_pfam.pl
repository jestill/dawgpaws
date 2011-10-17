#!/usr/bin/perl -w
#
#
# Download PFAM with wget
# prefix
my $url_src="http://pfam.sanger.ac.uk/family/hmm/";

@pfam_domains = (
    PF00075,
    PF00077,
    PF00078,
    PF00096,
    PF00385,
    PF00607,
    PF00552,
    PF00665,
    PF00872,
    PF01359,
    PF01498,
    PF01710,
    PF02022,
    PF02093,
    PF02813,
    PF02992,
    PF02994,
    PF03004,
    PF03017,
    PF03101,
    PF03108,
    PF03732,
    PF04195,
    PF04827,
    PF05699,
    PF06815,
    PF06817,
    PF07727,
    PF08284,
    PF09299,
    PF10536,
    PF10551,
    );

for my $pfam_domain (@pfam_domains) {
    my $wget_cmd = "wget ".$url_src.$pfam_domain;
    print STDERR $wget_cmd."\n";
    system ($wget_cmd);
    
    my $rename_cmd = "mv ".$pfam_domain." ".$pfam_domain.".hmm";
    print STDERR $rename_cmd."\n";
    system ($rename_cmd);
}
