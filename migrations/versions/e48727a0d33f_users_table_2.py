"""users table 2

Revision ID: e48727a0d33f
Revises: 
Create Date: 2022-12-08 16:58:42.362943

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = 'e48727a0d33f'
down_revision = None
branch_labels = None
depends_on = None


def upgrade():
    # ### commands auto generated by Alembic - please adjust! ###
    op.create_table('chemical',
    sa.Column('id', sa.Integer(), nullable=False),
    sa.Column('inchi', sa.String(), nullable=True),
    sa.PrimaryKeyConstraint('id')
    )
    # ### end Alembic commands ###


def downgrade():
    # ### commands auto generated by Alembic - please adjust! ###
    op.drop_table('chemical')
    # ### end Alembic commands ###
