#!/bin/bash

for x in a b c d e f g h i j k l m n o
do
  
echo "\$(OBJCLASS2$x): \$(CLASS2$x) "
echo "	@echo \"compling class file 2...\""
echo "	\$(CC) \$(CFLAG) -c \$(CLASS2$x) \$(INCLALL)"

done
