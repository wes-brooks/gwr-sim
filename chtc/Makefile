CHTCUTILS:=chtc-utils
export CHTCUTILS

submit-gaussian:
	./submit.sh

jobdirs-gaussian:
	find 'simdata-gaussian' -wholename 'simdata-gaussian/shared' -prune -o -not -name 'simdata-gaussian' | xargs rm -rf 
	rm -rf simoutput-gaussian
	rm -rf simdata-gaussian
	mkdir simdata-gaussian
	mkdir simdata-gaussian/shared
	cp -rf ../R-libs simdata-gaussian/shared/R-libs
	cp -rf ../data simdata-gaussian/shared/data
	cp -rf ../R/* simdata-gaussian/shared/
	cp -rf ../seeds.txt simdata-gaussian/shared
	./make-job-dirs.sh 1800 simrun-gaussian
	./mkdag \
	  --cmdtorun=beauty.r \
	  --data=simdata-gaussian \
	  --dagdir=simoutput-gaussian \
	  --pattern=result.txt \
	  --type=R \
	  --version=sl5-R-2.15.3
