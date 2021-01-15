#!/usr/bin/env python
import psycopg2

try:
   conn = psycopg2.connect("dbname='reg3sim' user='ajay' host='localhost' port='5432' password='adcirc'")
   cur = conn.cursor()

   cur.execute("""SELECT table_name FROM information_schema.tables WHERE table_schema = 'public'""")

   tables = cur.fetchall()

   f = open('/home/data/ingestProcessing/ajay/data/reg3sim_tables.txt','w')

   for table in tables:
       for row in table:
          f.write(row+'\n')

   f.close()

except (Exception, psycopg2.DatabaseError) as error:
   print(error)
finally:
   if cur is not None:
       cur.close()
   if conn is not None:
       conn.close()


