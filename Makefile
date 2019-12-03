
all: validate validate_with_json

WDL= $(shell find . -name '*.wdl')

all-tests := $(addsuffix .test, $(basename $(WDL)))
all-tests-json := $(addsuffix .json, $(all-tests))

validate: $(all-tests)
validate_with_json : $(all-tests) $(all-tests-json)

test : $(all-tests)

test_with_json: test $(all-tests-json)

%.test : %.wdl
	java -jar ~/womtool-46.jar validate $?

%.test.json : %.wdl.json
	java -jar /app/womtool-46.jar validate $(basename $? .json) -i $?