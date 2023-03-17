import yagmail
yag = yagmail.SMTP('sidescan.noreply@gmail.com', 'xzunrluyjhjzrpqr')
contents = ['This is the body, and here is just text',
            'You can find an audio file attached.']
yag.send('fred.s.kremer@gmail.com', 'xzunrluyjhjzrpqr')