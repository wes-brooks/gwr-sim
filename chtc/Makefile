CHTCUTILS:=chtc-utils
export CHTCUTILS


gaussian:
	find 'simdata-gaussian' -wholename 'simdata-gaussian/shared' -prune -o -not -name 'simdata-gaussian' | xargs rm -rf 
	rm -rf simoutput-gaussian
	rm -rf simdata-gaussian
	mkdir simdata-gaussian
	mkdir simdata-gaussian/shared
	cp -rf ../rlibs simdata-gaussian/shared/rlibs
	cp -rf ../R/* simdata-gaussian/shared/
	cp -rf ../seeds.txt simdata-gaussian/shared
	./make-job-dirs.sh gaussian 1800
	./mkdag \
	  --cmdtorun=simulation-gaussian.r \
	  --data=simdata-gaussian \
	  --dagdir=simoutput-gaussian \
	  --pattern=result.txt \
	  --type=R \
	  --version=sl5-R-2.15.3
	./submit.sh gaussian


binomial:
	find 'simdata-binomial' -wholename 'simdata-binomial/shared' -prune -o -not -name 'simdata-binomial' | xargs rm -rf 
	rm -rf simoutput-binomial
	rm -rf simdata-binomial
	mkdir simdata-binomial
	mkdir simdata-binomial/shared
	cp -rf ../R-libs simdata-binomial/shared/R-libs
	cp -rf ../R/* simdata-binomial/shared/
	cp -rf ../seeds.txt simdata-binomial/shared
	./make-job-dirs.sh binomial 900
	./mkdag \
	  --cmdtorun=simulation-binomial.r \
	  --data=simdata-binomial \
	  --dagdir=simoutput-binomial \
	  --pattern=result.txt \
	  --type=R \
	  --version=sl5-R-2.15.3
	./submit.sh binomial
	  


poisson:
	find 'simdata-poisson' -wholename 'simdata-poisson/shared' -prune -o -not -name 'simdata-poisson' | xargs rm -rf 
	rm -rf simoutput-poisson
	rm -rf simdata-poisson
	mkdir simdata-poisson
	mkdir simdata-poisson/shared
	cp -rf ../R-libs simdata-poisson/shared/R-libs
	cp -rf ../R/* simdata-poisson/shared/
	cp -rf ../seeds.txt simdata-poisson/shared
	./make-job-dirs.sh poisson 900
	./mkdag \
	  --cmdtorun=simulation-poisson.r \
	  --data=simdata-poisson \
	  --dagdir=simoutput-poisson \
	  --pattern=result.txt \
	  --type=R \
	  --version=sl5-R-2.15.3
	./submit.sh poisson
