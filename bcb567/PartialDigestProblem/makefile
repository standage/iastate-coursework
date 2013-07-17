JC = javac
CFLAGS = -cp .:mathCollection.jar
DRIVER = GenerateRestrictionMap
CLASSES = RestrictionMapper.class
TEST = RestrictionMapTest

proc:				$(TEST).class
				java $(CFLAGS) $(TEST)

$(TEST).class:			$(TEST).java
				$(JC) $(CFLAGS) $(TEST).java

oo:				$(DRIVER).class $(CLASSES)
				java $(CFLAGS) $(DRIVER)

$(DRIVER).class:		$(DRIVER).java $(CLASSES)
				$(JC) $(CFLAGS) $(DRIVER).java

RestrictionMapper.class:	RestrictionMapper.java
				$(JC) $(CFLAGS) RestrictionMapper.java

clean:
				rm *.class
