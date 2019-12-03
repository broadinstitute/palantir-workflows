
all: validate validate_with_json

WDL= $(shell find . -name '*.wdl')

VALIDATE=java -jar /app/womtool.jar validate

all-tests := $(addsuffix .test, $(basename $(WDL)))
all-tests-json := $(addsuffix .json, $(all-tests))

validate: $(all-tests)
validate_with_json : $(all-tests) $(all-tests-json)

test : $(all-tests)

test_with_json: test $(all-tests-json)

%.test : %.wdl
	$(VALIDATE)  $?

%.test.json : %.wdl.json
	$(VALIDATE) $(basename $? .json) -i $?