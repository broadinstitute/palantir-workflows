TEST_JSON= $(shell find test -name '*.json')

VALIDATE_WDL= $(shell find . -name '*.wdl' ! -path './test/*')

TEST=java -jar $(CROMWELL_JAR) run

VALIDATE=java -jar $(WOMTOOL_JAR) validate

all-tests := $(addsuffix .test, $(TEST_JSON))

all-validations := $(addsuffix .validate, $(VALIDATE_WDL))

.PHONY: all
all: test validate

.PHONY: test
test : $(all-tests)

.PHONY: validate
validate: $(all-validations)

.PHONY: %.wdl.validate
%.wdl.validate : 
	$(VALIDATE) $(basename $@)

.PHONY: %.json.test
%.json.test :
	$(TEST) $(subst _json/,.wdl, $(dir $(basename $@))) -i $(basename $@)
