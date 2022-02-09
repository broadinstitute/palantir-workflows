SHELL=/bin/bash -o pipefail

TEST_JSON= $(shell find test -name '*.json')

VALIDATE_WDL= $(shell find . -name '*.wdl' ! -path './test/*')

TEST=java -jar $(CROMWELL_JAR) run

VALIDATE=java -jar $(WOMTOOL_JAR) validate

all-tests := $(addsuffix .test, $(TEST_JSON))

all-validations := $(addsuffix .validate, $(VALIDATE_WDL))

all-validations-with-json := $(addsuffix .validate, $(TEST_JSON))

.PHONY: all
all: test validate

.PHONY: test
test: $(all-tests)

.PHONY: validate
validate: $(all-validations)

.PHONY: validate-with-json
validate-with-json: $(all-validations-with-json)

.PHONY: %.wdl.validate
%.wdl.validate: 
	$(VALIDATE) $(basename $@)

.PHONY: %.json.test
%.json.test:
	mkdir -p logs
	$(TEST) $(subst _json/,.wdl, $(dir $(basename $@))) -i $(basename $@) -o test_options.json 2>&1 | tee logs/$(notdir $(subst _json/,.wdl, $(dir $(basename $@))))_with_$(notdir $(basename $@)).log ; \
	test_result=$${PIPESTATUS[0]} ; \
	if [ $${test_result} -ne 0 ] ; then \
	    echo "=================================== FAILURE ===================================" ; \
	    echo "Test failed for $(notdir $(subst _json/,.wdl, $(dir $(basename $@)))) with $(notdir $(basename $@)). Check the output logs (uploaded as artifacts) for more information." ; \
	    echo "=================================== FAILURE ===================================" ; \
	fi ; \
	exit $${test_result}

.PHONY: %.json.validate
%.json.validate:
	$(VALIDATE) $(subst _json/,.wdl, $(dir $(basename $@))) -i $(basename $@)
