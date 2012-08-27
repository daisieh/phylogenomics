my $pid = open (READ, "-|");
defined($pid)           || die "can't fork: $!";
if ($pid) {             # parent
	print "parent here...\n";
	while (my $line = <READ>) {
		my $x = select (STDERR);
		$| = 1;
		select ($x);
		print STDERR "yep $line";
						# do something interesting
	}
	close(READ)  || warn "kid exited $?";
} else {                # child
	print "kid here...\n";
	($EUID, $EGID) = ($UID, $GID); # suid only
 	my $pid2 = open (KID, "HYPHYSP BASEPATH= /Users/daisie/Documents/Work/Sandbox/bio_scripts/test.bf |") || die "can't exec program: $!";
# 	open (OUTFILE, ">~/Desktop/outfile");
#	pipe (KID, READ);
# 	my $x = select (OUTFILE);
# 	$| = 1;
# 	select ($x);
	print READ "hola\n";
	waitpid $pid2, 0;
# 	close OUTFILE;
	#exec($program, @options, @args) || die "can't exec program: $!";
	# NOTREACHED
}
