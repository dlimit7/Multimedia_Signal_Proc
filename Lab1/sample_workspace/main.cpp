#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
	if (argc > 1)
	{
		char str[BUFSIZ];
		sscanf(argv[1], "%s", str);
		printf("Hello, my name is %s\n", str);
	}
 
  return 0;
}


